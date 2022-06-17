# requires bedtools to be installed on the command line

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  mergeOrder <- c(here('extractedData', 'HINT', 'consensusPeakFiles', 'RAhigh.consensusSummits.bed'),  
                  here('extractedData', 'HINT', 'consensusPeakFiles', 'TGFbhigh.consensusSummits.bed'),  
                  here('extractedData', 'HINT', 'consensusPeakFiles', 'EtOH.consensusSummits.bed'))  
  inputConsensusPeaksFolder <-  here('extractedData', 'HINT', 'consensusPeakFiles')
  outputMergedPeaksFile <- sprintf('"%s"', here('extractedData', 'HINT', 'mergedConsensusPeakFiles', 'allCondsMergedSummitWindows.bed'))
} else {
  inputConsensusPeaksFolder <- cmdargs[1]
  outputMergedPeaksFile <- sprintf('"%s"', cmdargs[2])
  
  mergeOrder <- c(paste0(inputConsensusPeaksFolder, '/', 'RAhigh.consensusSummits.bed'),  
                  paste0(inputConsensusPeaksFolder, '/', 'TGFbhigh.consensusSummits.bed'),  
                  paste0(inputConsensusPeaksFolder, '/', 'EtOH.consensusSummits.bed'))
}



mergeOrder <- map(mergeOrder, function (x) sprintf('"%s"', x)) #sprintf functions needed for command line parsing


tempBedFile.forMerging  <- sprintf('"%s"', paste0(inputConsensusPeaksFolder, '/', 'tempMergedFile.bed'))
tempBedFile.forMerging2 <- sprintf('"%s"', paste0(inputConsensusPeaksFolder, '/', 'tempMergedFile2.bed'))
tempBedFile.forDiffs    <- sprintf('"%s"', paste0(inputConsensusPeaksFolder, '/', 'tempDiffPeaks.bed'))
tempBedFile.forSorting  <- sprintf('"%s"', paste0(inputConsensusPeaksFolder, '/', 'tempMergedSortedFile.bed'))
#step 1: merge overlapping segments in the first BED file
system(paste('bedtools', 'merge', '-c', '4,5', '-o', 'collapse,max', '-i', mergeOrder[1], '>', tempBedFile.forSorting))  
#next steps: in order, add new conditions to existing peaks, only if they don't already overlap a current peak
for (fileToMergeIn in mergeOrder[2:length(mergeOrder)]) {
  system(paste('bedtools', 'merge', '-c', '4,5', '-o', 'collapse,max', '-i', fileToMergeIn, '>', tempBedFile.forMerging2))
  system(paste('bedtools', 'intersect', '-v', '-sorted', '-a', tempBedFile.forMerging2, '-b', tempBedFile.forSorting, '>', tempBedFile.forDiffs))
  system(paste('cat', tempBedFile.forDiffs, tempBedFile.forSorting, ">", tempBedFile.forMerging))
  system(paste('sort', '-k1,1', '-k2,2n', tempBedFile.forMerging, '>', tempBedFile.forSorting))
}
system(paste('cp', tempBedFile.forSorting, outputMergedPeaksFile))
system(paste('rm', tempBedFile.forMerging, tempBedFile.forMerging2, tempBedFile.forDiffs, tempBedFile.forSorting))
