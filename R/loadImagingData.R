
#' Loads the images for a set of slide regions
#'
#' Just a convenience function that gets directories to all the imaging data on the server.
#' @param slideRegionList A character vector of slide_regions to load in. Slide and region should be separated by an underscore.
#' @param markers qcname .csv file with qc variables
#' @param adjustedImageDir directory of the AF adjusted images.
#' @param maskdir directory of where the mask files are stored.
#' @return Returns data frame with image directories.
#' @export
#' @importFrom utils read.csv
loadImagingData = function(slideRegionList, markers, adjustedImageDir='/media/disk2/atlas_mxif/BALANCED', maskdir='/media/disk2/atlas_mxif/BALANCED'){
  # This gets regions whose images exist for each slide
  #region = unlist(lapply(slideList, function(slide) gsub("^DAPI_ADJ_|\\.jpg", "", list.files(file.path(adjustedImageDir, slide), pattern='DAPI_ADJ_*'))) )
  
  
  #idVars = c('Slide', 'Region', 'slideRegion')
  #sr = threshold[!duplicated(threshold$slideRegion), idVars]
  sr = data.frame(slide_region = slideRegionList, 'region'=gsub('.*_region_', 'region_', slideRegionList), 'slide'=gsub('_region.*', '', slideRegionList))
  sr[, markers ] = sapply(markers, function(marker) file.path(adjustedImageDir, sr$slide, sprintf('%s_ADJ_%s.jpg', marker, sr$region) ) )
  #sr[,'mask'] = file.path('/media/disk2/atlas_mxif/T_cell_quantification/Images/Adjusted Images', sr$Slide, sprintf('%s_%02d_TISSUE_MASK.tif', sr$Slide, as.numeric(gsub('region_', '', sr$region)) ) )
  sr[,'tumorMask'] = file.path(maskdir, sr$slide, sprintf('Tumor_mask_%s.png', sr$region) )
  #sr[,'epiMask'] = file.path('/media/disk2/atlas_mxif/T_cell_quantification/Images/', sr$Slide, sprintf('%s_%s_epi_mask.png', sr$Slide, sr$region) )
  #sr[,'strMask'] = file.path('/media/disk2/atlas_mxif/T_cell_quantification/Images/', sr$Slide, sprintf('%s_%s_stroma_mask.png', sr$Slide, sr$region) )
  
  sr[,markers][!apply(sr[,markers], 2, file.exists)] = NA
  sr
}
