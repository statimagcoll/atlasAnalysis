
#' Hierarchical clustering using naive Bayes
#'
#' Old script that combines separate clustering algorithms together in a hierarchical approach. Should be revamped before use.
#' @param slidedata data frame of the markers for each cell.
#' @param markers data frame indicating the hierarchical structure of the cluster types.
#' @param probthr Threshold to define marker positive cells
#' @importFrom combinat permn
#' @importFrom mclust Mclust priorControl
#' @importFrom utils capture.output
#' @importFrom stats predict sd
HCnaive = function(slidedata, markers, probthr=0.8){
clusters = grep('cl[0-9]$', names(markers), value=TRUE)
for(cl in clusters){
  cllevel = as.numeric(gsub('cl', '', cl))
  classes = unique(markers[,cl])
  clResults = list()
  for(clas in classes){
    if(cllevel>1){
      upperclas = markers[which(markers[,cl] %in% clas)[1], paste0('cl', cllevel-1)]
      pupperclas = paste(paste0('cl', cllevel-1), upperclas, sep='_')
    }
    markernames = paste0('Median_Cell_', unlist(markers[which(markers[,cl] %in% clas)[1], paste0(cl, 'markers') ]))
    markervalues = unlist(markers[which(markers[,cl] %in% clas)[1], paste0(cl, 'markerValue')])
    result = apply(slidedata[,c(markernames), drop=FALSE], 2, function(x, modelSubset){
      n = length(x)
      nzero = sum(x<=0)
      zeros = (x<=0)
      # only modelSubset are used to build the model
      xsubset = x[which(!zeros & modelSubset)]
      invisible(capture.output(mixtmod <- tryCatch(mclust::Mclust(xsubset, G=2, model='E', prior=priorControl()), error=function(e){NA})))
      ans = rep(NA, n)
      if(!is.na(mixtmod[1])){
        maxind = which.max(mixtmod$parameters$mean)
        # all things that are zero get probability zero of being in the higher mixture group
        ans[which(zeros)] = 0
        ans[which(!zeros)] = predict(mixtmod, newdata = x[which(!zeros)])$z[,maxind]
      }
      ans
    }, modelSubset=if(cllevel>1) ans[,pupperclas]>probthr else (rep(TRUE, nrow(slidedata) ) ) )
    # if markervalues ==2 return probability of being in the higher mean cluster, else return probability of being in the lower mean cluster.
    # label is mean of log(probabilities)
    missings = apply(is.na(result), 2, all)
    if(any(missings)){
      warning('Clustering failed for marker ', markernames[missings])
    }
    clResults[[clas]] = apply(sapply(1:length(markernames), function(ind){ if(all(is.na(result[,ind]))){ ans = NULL} else {ans = if(markervalues[ind]==2) result[,ind] else 1-result[,ind]; log(ans); }; ans}, simplify='array'), 1, mean)
  } # end classes

  # convert back to probability scale
  probs = exp(do.call(cbind, clResults))
  # probabilities with numerically zero values create issues
  probs[probs==0 ] = 10^{-8}
  if(cllevel>1){
    upperclasses = markers[,paste0('cl', cllevel-1)]
    classesrows = markers[,cl]
    # scale probabilities to sum to one within each upper cluster
    probs = lapply(unique(upperclasses), function(uclass){columns =  which(colnames(probs) %in% classesrows[ which(upperclasses %in% uclass)]);  sweep(probs[,columns], 1, rowSums(probs[,columns]), FUN='/')} )
    probs = do.call(cbind, probs)
    colnames(probs) = paste(cl, classes, sep='_')
    ans[, colnames(probs)] = probs
  } else {
    # scale probabilities to sum to one
    ans = as.data.frame(sweep(probs, 1, rowSums(probs), FUN='/'))
    names(ans) = paste(cl, classes, sep='_')
  }
} # end clusters
# return input data and computed probabilities
list(data=slidedata, z=ans)
}






### This function is not complete, need to write it for level 2 clustering
HCmulti = function(slidedata, markers, probthr=0.5){
  clusters = grep('cl[0-9]$', names(markers), value=TRUE)
  for(cl in clusters){
    cllevel = as.numeric(gsub('cl', '', cl))
      if(cllevel>1){
        upperclasses = unique(markers[, paste0('cl', cllevel-1)])
        for(uclas in upperclasses){
          classes = unique(markers[markers[, paste0('cl', cllevel-1)]==uclas,cl])
          nclasses = length(classes)
          inds = ans[[paste0('cl', cllevel-1)]]$z[,uclas]>probthr
          markernames = unique(unlist(markers[ markers[, paste0('cl', cllevel-1)]==uclas,paste0(cl, 'markers')]))

          zeromarkers = markernames[ apply(slidedata[,markernames]==0, 2, all)]
          mclmod = Mclust(slidedata[inds,markernames[! markernames %in% zeromarkers]], G=nclasses, prior=priorControl(), modelNames='EEI')
          # subtracts the minimum mean for that marker across clusters and scales by the SD of the marker
          scmeans = t(scale(t(mclmod$parameters$mean), center=apply(mclmod$parameters$mean, 1, min), scale=apply(slidedata[inds,markernames[! markernames %in% zeromarkers]], 2, sd)))
          scmeans = sweep(scmeans, 1, apply(scmeans, 1, max), '/')
          scmeans = rbind(scmeans, matrix(0, ncol=ncol(scmeans), nrow=length(zeromarkers)) )
          rownames(scmeans) = c(rownames(mclmod$parameters$mean), zeromarkers)
          markervals = as.matrix(scmeans); markervals[,] = NA; colnames(markervals) = classes
          for(clas in classes){
            clmarkers = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markers')])
            clmarkervalues = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markerValue')])
            markervals[clmarkers, clas] = clmarkervalues
          }
          combs = as.data.frame(do.call(rbind, combinat::permn(1:nclasses)))
          combs$dist = apply(combs, 1, function(inds){dist = sum(colMeans(abs(markervals[,inds] - scmeans), na.rm=TRUE)); dist})
          # add column names
          colnames(mclmod$parameters$mean) = colnames(mclmod$z) = classes[ unlist(combs[which.min(combs$dist), 1:nclasses])]
          # add indices
          mclmod$inds = inds
          mclmod$scmean = scmeans
          mclmod$dist = min(combs$dist)
          ans = c(ans, list(mclmod))
        }
      } else {
      classes = unique(markers[,cl])
      nclasses = length(classes)
      markernames = unique(unlist(markers[,paste0(cl, 'markers')]))

      zeromarkers = markernames[ apply(slidedata[,markernames]==0, 2, all)]
      mclmod = Mclust(slidedata[,markernames[! markernames %in% zeromarkers]], G=nclasses, prior=priorControl(), modelNames='EEI')
      # subtracts the minimum mean for that marker across clusters and scales by the SD of the marker
      scmeans = t(scale(t(mclmod$parameters$mean), center=apply(mclmod$parameters$mean, 1, min), scale=apply(slidedata[,markernames[! markernames %in% zeromarkers]], 2, sd)))
      scmeans = sweep(scmeans, 1, apply(scmeans, 1, max), '/')
      scmeans = rbind(scmeans, matrix(0, ncol=ncol(scmeans), nrow=length(zeromarkers)) )
      rownames(scmeans) = c(rownames(mclmod$parameters$mean), zeromarkers)
      markervals = as.matrix(scmeans); markervals[,] = NA; colnames(markervals) = classes
      for(clas in classes){
        clmarkers = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markers')])
        clmarkervalues = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markerValue')])
        markervals[clmarkers, clas] = clmarkervalues
      }
      combs = as.data.frame(do.call(rbind, combinat::permn(1:nclasses)))
      combs$dist = apply(combs, 1, function(inds2){dist = sum(colMeans(abs(markervals[,inds2] - scmeans), na.rm=TRUE)); dist})
      # add column names
      colnames(mclmod$parameters$mean) = colnames(mclmod$z) = classes[ unlist(combs[which.min(combs$dist), 1:nclasses])]
      mclmod$scmean = scmeans
      mclmod$dist = min(combs$dist)
      ans = list('cl1'=mclmod)
      }
  } # end clusters
  # return input data and computed probabilities
  ans
}

# markers should be a subset of the markers data frame of interest
assignClasses = function(mclmod, markers, cl='cl1'){
  classes = unique(markers[,cl])
  nclasses = length(classes)
  markernames = unique(unlist(markers[,paste0(cl, 'markers')]))
  zeromarkers = markernames[ ! markernames %in% rownames(mclmod$parameters$mean)]

  #scmeans = t(scale(t(mclmod$parameters$mean), center=apply(mclmod$parameters$mean, 1, min), scale=apply(mclmod$data, 2, sd)))
  #scmeans = sweep(scmeans, 1, apply(scmeans, 1, max), '/')
  scmeans = t(scale(t(mclmod$parameters$mean), center=apply(mclmod$parameters$mean, 1, min), scale=apply(mclmod$parameters$mean, 1, max) ) )
  scmeans = rbind(scmeans, matrix(0, ncol=ncol(scmeans), nrow=length(zeromarkers)) )
  rownames(scmeans) = c(rownames(mclmod$parameters$mean), zeromarkers)
  markervals = as.matrix(scmeans); markervals[,] = NA; colnames(markervals) = classes
  for(clas in classes){
    clmarkers = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markers')])
    clmarkervalues = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markerValue')])
    markervals[clmarkers, clas] = clmarkervalues
  }
  combs = as.data.frame(do.call(rbind, combinat::permn(1:nclasses)))
  combs$dist = apply(combs, 1, function(inds){dist = sum(colMeans(abs(markervals[,inds] - scmeans), na.rm=TRUE)); dist})
  # add column names
  colnames(mclmod$parameters$mean) = colnames(mclmod$z) = classes[ unlist(combs[which.min(combs$dist), 1:nclasses])]
  # add indices
  mclmod$scmean = scmeans
  mclmod$dist = min(combs$dist)
  mclmod
}

# markers should be a subset of the markers data frame of interest
assignClasses.nmf = function(basis, coefs, markers, cl='cl1'){
  classes = unique(markers[,cl])
  nclasses = length(classes)
  markernames = unique(unlist(markers[,paste0(cl, 'markers')]))
  zeromarkers = markernames[ ! markernames %in% colnames(coefs)]

  #scmeans = t(scale(t(mclmod$parameters$mean), center=apply(mclmod$parameters$mean, 1, min), scale=apply(mclmod$data, 2, sd)))
  #scmeans = sweep(scmeans, 1, apply(scmeans, 1, max), '/')
  scmeans = t(scale(coefs, center=apply(coefs, 2, min), scale=apply(coefs, 2, max) ) )
  scmeans = rbind(scmeans, matrix(0, ncol=ncol(scmeans), nrow=length(zeromarkers)) )
  rownames(scmeans) = c(colnames(coefs), zeromarkers)
  markervals = as.matrix(scmeans); markervals[,] = NA; colnames(markervals) = classes
  for(clas in classes){
    clmarkers = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markers')])
    clmarkervalues = unlist(markers[ which(markers[,cl] %in% clas)[1], paste0(cl, 'markerValue')])
    markervals[clmarkers, clas] = clmarkervalues
  }
  combs = as.data.frame(do.call(rbind, combinat::permn(1:nclasses)))
  combs$dist = apply(combs, 1, function(inds){dist = sum(colMeans(abs(markervals[,inds] - scmeans), na.rm=TRUE)); dist})
  # add column names
  rownames(coefs) = colnames(basis) = classes[ unlist(combs[which.min(combs$dist), 1:nclasses])]
  basis
}
