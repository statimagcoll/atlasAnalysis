
#' Loads the quantifed cell level PCA atlas data into memory
#'
#' Just a convenience function so that we're all using the same data set
#' @param dataname A data frame with the cell level data
#' @param qcname .csv file with qc variables
#' @return Returns data frame with the data merged.
#' @export
#' @importFrom utils read.csv
loadQuantifiedPCA = function(dataname='/media/disk2/atlas_mxif/data/colon map batch 1-3 position corrected w immune labels.rds', qcname='/media/disk2/atlas_mxif/data/Colon MAP QC.csv'){
  sc = readRDS(dataname)
  tts = read.csv(qcname)

  # Eliot's data set
  tts$SlideID = sc$SlideID[match(tts$TissueID, sc$TissueID)]
  sc$Slide_Region = paste(sc$SlideID, sc$region, sep='_')
  tts$sampled = tts$SlideID %in% unique(sc$SlideID)
  sc = merge(sc, tts, all.x=TRUE)
  sc
}
