library("tidyverse")
library(chromVAR)
library(GenomicRanges)
library(here)


cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peak.set.to.use <- "all" # should be one of [all, differential, or an absolute path to a bed file]
  peakWidth <- 150 # we used a fixed window size of 150bp (~1 nucleosome), centered at MACS2 summits.
  # Since they're non-overlapping, we increase the size to 225 for motif-analysis
  # to reduce the chance of missing motifs in the boundary regions between consecutive peaks
  center.peaks <- F
} else {
  peak.set.to.use <- cmdargs[1] # should be one of [all, differential, or an absolute path to a bed file]
  peakWidth       <- as.numeric(cmdargs[2])
  outputFile      <- cmdargs[3]
  center.peaks <- !(cmdargs[4] == "VariableWidth")
  print(center.peaks)
  print(outputFile)
}

sampleMetadata  <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'il1a_sampleMetadata.txt'), col_names=TRUE)
tn5insertionPointsBed <- paste(sampleMetadata$ATAC_analysisDir, "/", sampleMetadata$sampleID, '.Tn5_insertion_points.tagAlign.gz', sep="")
sampleMetadata  <- mutate(sampleMetadata, tn5insertionPointsBed=tn5insertionPointsBed)
relevantMetadata <- sampleMetadata




if (peak.set.to.use == "all") {
  mergedPeaksFile <- here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'mergedConsensusPeakFiles', 'allCondsMergedSummitWindows_mergedist50.bed')
  outputFile <- here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'atacFragmentCountsAllCondsMergedSummitWindows.rds')
} else if (peak.set.to.use == "differential") {
  mergedPeaksFile <- here('extractedData', 'differentialAtacPeaks_forIGV.bed')
  outputFile <- here('extractedData', 'atacFragmentCountsAllCondsDifferentialPeaks.rds')
} else if (peak.set.to.use == "differential_merged") {
  mergedPeaksFile <- here('extractedData', 'differentialAtacPeaks_merged_forIGV.bed')
  outputFile <- here('extractedData', 'atacFragmentCountsAllCondsMergedDifferentialPeaks.rds')
  center.peaks <- F
} else {
  mergedPeaksFile <- peak.set.to.use
}


peakSetGRanges <- getPeaks(mergedPeaksFile)
if (center.peaks) {
  peakSetGRanges <- resize(peakSetGRanges, width = peakWidth, fix = "center")
}
conditionsToCount <- unique(relevantMetadata$condition)
tn5bedFileList <- filter(relevantMetadata, condition %in% conditionsToCount)$tn5insertionPointsBed
sampleIDs <- filter(relevantMetadata, condition %in% conditionsToCount)$sampleID
# Caution: do not use bams as input! for bams, chromVar assigns the entire range covered by a read pair to the peaks it overlaps
fragmentCounts <- getCounts(tn5bedFileList, 
                            peakSetGRanges, 
                            paired =  FALSE, 
                            by_rg = FALSE, 
                            format = "bed", 
                            colData = DataFrame(celltype = sampleIDs))

saveRDS(fragmentCounts, file=outputFile)
