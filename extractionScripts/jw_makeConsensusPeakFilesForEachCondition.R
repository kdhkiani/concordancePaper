#########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  sampleMetadata <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'il1a_sampleMetadata.txt'), col_names=TRUE)
  outputFolder       <- here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'consensusPeakFiles')
} else {
  sampleMetadataFile <- cmdargs[1]
  outputFolder       <- cmdargs[2]
}

# In future work, can potentially add the ability to map the sample metadata to cloud-based storage of the ATAC pipeline output files below
sampleMetadata <- mutate(sampleMetadata, macs2summits = paste(ATAC_analysisDir, '/macs2_output_windowSize150/', `sampleID`, '_summits.blacklistFiltered.bed', sep=""))
conditionsToFindConsensusPeaksFor <- sampleMetadata %>% dplyr::select(condition) %>% unique() %>% as_vector()
peakMergeWindowSize <- 150
#nReps <- 3
summitWindowExtensionRadius <- 75

# function for determining consensus peak set for a condition with three replicates
calculateAndWriteConsensusPeakSet <- function(conditionName, sampleMetadata, peakMergeWindowSize, summitWindowExtensionRadius, nReps) {
  
  conditionMetadataAllReps <- filter(sampleMetadata, condition==conditionName)
  
  summit_colnames <- c('chr', 'startPos', 'endPos', 'peakID', 'negLogPval')
  summits1 <- read_tsv(conditionMetadataAllReps$macs2summits[1], col_names = summit_colnames)
  
  summitsAll <- summits1
  summitsAll <- separate(summitsAll, peakID, c('sampleName', 'peakNum'), sep="_peak_")
  summitsAllSorted <- arrange(summitsAll, chr, startPos)
  
 
  
  selectedSummitsSorted <- summitsAllSorted
  # For final summit file, we want each 
  selectedSummitsSorted$startPos <- selectedSummitsSorted$startPos - summitWindowExtensionRadius
  selectedSummitsSorted$endPos   <- selectedSummitsSorted$endPos + summitWindowExtensionRadius
  
  selectedSummitsToWrite <- unite(selectedSummitsSorted, 'peakID', c('sampleName', 'peakNum'), sep = "_peak_", remove = TRUE)
  
  outputFile <- paste0(outputFolder, '/', gsub(' ', '-', conditionName), '.consensusSummits.bed')
  write_tsv(selectedSummitsToWrite, outputFile, col_names=FALSE)
}

# Helper function for calculateAndWriteConsensusPeakSet. This function requires SORTED input location vectors.
getIndicesWithinWindowSize <- function(summitStartLocs, currentIndex, windowSize, nTotalPeaks) {
  currentSummitLoc <- summitStartLocs[currentIndex]
  indicesWithinWindow <- c(currentIndex)
  nextIndex <- currentIndex + 1
  nextIndexWithinWindow <- TRUE
  while (nextIndexWithinWindow & nextIndex <= nTotalPeaks) {
    nextSummitLoc <- summitStartLocs[nextIndex]
    nextIndexWithinWindow <- abs(nextSummitLoc - currentSummitLoc) <= windowSize
    if (nextIndexWithinWindow) {
      indicesWithinWindow <- c(indicesWithinWindow, nextIndex)
    }
    nextIndex <- nextIndex + 1
  }
  return(indicesWithinWindow)
}

# Script runs in loop below. 
for (condition in conditionsToFindConsensusPeaksFor) {
  calculateAndWriteConsensusPeakSet(condition, sampleMetadata, peakMergeWindowSize, summitWindowExtensionRadius, nReps)
}





