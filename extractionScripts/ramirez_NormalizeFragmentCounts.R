############################## user-defined parameters go in this block of code ##################################
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCounts <- readRDS(here('extractedData', 'Ramirez', 'ATACseq', 'atacFragmentCountsAllCondsMergedSummitWindows.rds'))
  # ran
  mergedPeaksFile <- here('extractedData', 'Ramirez', 'ATACseq', 'mergedConsensusPeakFiles', 'allCondsMergedSummitWindows_mergedist50.bed')
  outputFile <- here('extractedData', 'Ramirez', 'ATACseq', 'ramirez_AtacPeaksMergedSummitWindows_mergedist50.tsv')
  outputBed <- here('extractedData', 'Ramirez', 'ATACseq', 'ramirez_AtacPeaksMergedSummitWindows_mergedist50.bed')
} else {
  fragmentCounts <- readRDS(cmdargs[1])  
  outputFile <- cmdargs[2]
}
peaks <- read_tsv(mergedPeaksFile, col_names = F)
peaks <- peaks %>% 
  dplyr::mutate(peakName = paste(peaks$X1,peaks$X2,peaks$X3,sep = ":"))

sampleMetadata        <- read_tsv(here('extractedData', 'Ramirez', 'hl60_sampleMetadata.txt'), col_names=TRUE)
tn5insertionPointsBed <- paste(sampleMetadata$ATAC_analysisDir, "/", sampleMetadata$sampleID, '.Tn5_insertion_points.tagAlign.gz', sep="")
sampleMetadata        <- mutate(sampleMetadata, tn5insertionPointsBed=tn5insertionPointsBed)
filtMetadataOrdered   <- sampleMetadata
conditionsToCount     <- unique(filtMetadataOrdered$condition)
nConditions           <- length(conditionsToCount)
nReplicates           <- 3

countMatrix           <- as.matrix(counts(fragmentCounts))
colnames(countMatrix) <- paste(filtMetadataOrdered$sampleID, filtMetadataOrdered$condition, sep = "_")
nSamples              <- length(colnames(countMatrix))
nPeaks                <- nrow(countMatrix)

sampleSizeFactors     <- nSamples * colSums(countMatrix) / sum(countMatrix)
normalizedCountMatrix <- t(t(countMatrix) / sampleSizeFactors)  # normalization is based on the number of reads in peaks

out_tib <- as_tibble(normalizedCountMatrix) 

out_tib <- out_tib %>% 
  dplyr::mutate(peak = peaks$peakName) %>% 
  dplyr::relocate(peak)

write_tsv(x = out_tib, file = outputFile)

out_bed <- out_tib %>% 
  dplyr::select(peak) %>% 
  tidyr::separate(peak, into = c("chrom", "start", "end"), sep = ":")

write_tsv(x = out_bed, file = outputBed, col_names = FALSE)
