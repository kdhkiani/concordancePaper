#########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  inputDir           <- here('extractedData', 'HINT')
  outputFolder       <- here('extractedData', 'HINT', 'consensusPeakFiles', '')
} else {
  inputDir           <- cmdargs[1]
  outputFolder       <- cmdargs[2]
}

# In future work, can potentially add the ability to map the sample metadata to cloud-based storage of the ATAC pipeline output files below
conditionsToFindConsensusPeaksFor <- c("EtOH", "RAhigh", "TGFbhigh")
peakMergeWindowSize <- 150
nReps <- 3

# function for determining consensus peak set for a condition with three replicates
calculateAndWriteConsensusPeakSet <- function(conditionName, peakMergeWindowSize, nReps) {
  stopifnot(nReps == 3)
  
  conditionMetadataAllReps <- list.files(path = inputDir, pattern = conditionName)
  
  summit_colnames <- c('chr', 'startPos', 'endPos', 'peakID', 'val', 'skip', 'na')
  summits1 <- read_tsv(here(inputDir, conditionMetadataAllReps[1]), col_names = summit_colnames)
  summits2 <- read_tsv(here(inputDir, conditionMetadataAllReps[2]), col_names = summit_colnames)
  summits3 <- read_tsv(here(inputDir, conditionMetadataAllReps[3]), col_names = summit_colnames)
  
  summits1 <- summits1 %>% dplyr::mutate(peakID = paste0(conditionMetadataAllReps[1],"_peak_", peakID))
  summits2 <- summits2 %>% dplyr::mutate(peakID = paste0(conditionMetadataAllReps[2],"_peak_", peakID))
  summits3 <- summits3 %>% dplyr::mutate(peakID = paste0(conditionMetadataAllReps[3],"_peak_", peakID))
  
  summitsAll <- rbind(summits1, summits2, summits3)
  summitsAll <- separate(summitsAll, peakID, c('sampleName', 'peakNum'), sep="_peak_")
  summitsAllSorted <- arrange(summitsAll, chr, startPos)
  
  summitsAllSorted <- summitsAllSorted %>% dplyr::mutate(summitLocs = startPos + round(((endPos - startPos)/2),0))
  
  # In this step, we iterate over all peaks to identify trios within 150 bp. They must include a peak from each of the 3 samples to count as a trio.
  nTotalPeaks <- nrow(summitsAll)
  # eligiblePeaks have not been marked as part of trios
  eligiblePeaks <- rep(TRUE, nTotalPeaks)
  # trioCenters are the centers of trios of peaks
  trioCenters <- rep(FALSE, nTotalPeaks)
  summitLocs <- summitsAllSorted$summitLocs
  summitSampleNames <- summitsAllSorted$sampleName
  for (ii in 1:nTotalPeaks) {
    if (eligiblePeaks[ii]) {
      indicesWithinWindow <- getIndicesWithinWindowSize(summitLocs, ii, peakMergeWindowSize, nTotalPeaks)
      if (length(unique(summitSampleNames[indicesWithinWindow])) == nReps) {
        includedSampleNames <- c()
        includedTrioIndices <- c()
        for (jj in indicesWithinWindow) {
          if (!(summitSampleNames[jj] %in% includedSampleNames)) {
            eligiblePeaks[jj] <- FALSE
            includedTrioIndices <- c(includedTrioIndices, jj)
            includedSampleNames <- c(includedSampleNames, summitSampleNames[jj])
          }
        }
        trioCenters[includedTrioIndices[2]] <- TRUE
      }
    }
  }
  selectedTrioSummits <- summitsAllSorted[trioCenters, ]
  
  # In this step, we iterate over all remaining non-Trio peaks to identify duos within 150 bp. They must include a peak from two different samples to count as a duo.
  # The peak with the stronger p-value is kept for duos.
  summitsAllSortedTriosRemoved <- summitsAllSorted[eligiblePeaks,]
  nPeaksForDuos <-  nrow(summitsAllSortedTriosRemoved)
  # eligiblePeaks have not been marked as part of duos
  eligiblePeaks <- rep(TRUE, nPeaksForDuos)
  # selectedDuoMembers are the member of a duo that has a stronger p value
  selectedDuoMembers <- rep(FALSE, nPeaksForDuos)
  summitLocs <- summitsAllSortedTriosRemoved$startPos
  summitSampleNames <- summitsAllSortedTriosRemoved$sampleName
  summitNegLogPvals <- summitsAllSortedTriosRemoved$val
  for (ii in 1:nPeaksForDuos) {
    if (eligiblePeaks[ii]) {
      indicesWithinWindow <- getIndicesWithinWindowSize(summitLocs, ii, peakMergeWindowSize, nPeaksForDuos)
      if (length(unique(summitSampleNames[indicesWithinWindow])) >= 2) {
        includedSampleNames <- c()
        includedDuoIndices <- c()
        for (jj in indicesWithinWindow) {
          if (!(summitSampleNames[jj] %in% includedSampleNames)) {
            eligiblePeaks[jj] <- FALSE
            includedDuoIndices <- c(includedDuoIndices, jj)
            includedSampleNames <- c(includedSampleNames, summitSampleNames[jj])
          }
        }
        # select duo member with higher negLog p value
        if (summitNegLogPvals[includedDuoIndices[1]] >= summitNegLogPvals[includedDuoIndices[2]]) {
          selectedDuoMembers[includedDuoIndices[1]] <- TRUE
        } else {
          selectedDuoMembers[includedDuoIndices[2]] <- TRUE
        }
      }
    }
  }
  selectedDuoSummits <- summitsAllSortedTriosRemoved[selectedDuoMembers, ]
  
  selectedSummits <- rbind(selectedDuoSummits, selectedTrioSummits)
  selectedSummitsSorted <- arrange(selectedSummits, chr, startPos)
  # For final summit file, we want each 
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
  calculateAndWriteConsensusPeakSet(condition, peakMergeWindowSize, nReps)
}
