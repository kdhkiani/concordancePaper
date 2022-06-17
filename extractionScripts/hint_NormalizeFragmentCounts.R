############################## user-defined parameters go in this block of code ##################################
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCounts <- readRDS(here('extractedData', 'HINT', 'hint_atacFragmentCountsAllCondsMergedSummitWindows.rds'))
  # ran
  mergedPeaksFile <- here('extractedData', 'HINT', 'mergedConsensusPeakFiles', 'allCondsMergedSummitWindows.bed')
  outputFile <- here('extractedData', 'HINT', 'hint_AtacPeaksMergedSummitWindows.tsv')
  outputBed <- here('extractedData', 'HINT', 'hint_AtacPeaksMergedSummitWindows.bed')
} else {
  fragmentCounts <- readRDS(cmdargs[1])  
  outputFile <- cmdargs[2]
}
peaks <- read_tsv(mergedPeaksFile, col_names = F)
peaks <- peaks %>% 
  dplyr::mutate(peakName = paste(peaks$X1,peaks$X2,peaks$X3,sep = ":"))

sampleMetadata  <- read_tsv(here('sampleMetadata_SI2-SI4.txt'), col_names=TRUE) #column of file location altered to reflect directory structure in removable HD
tn5insertionPointsBed <- paste(sampleMetadata$ATAC_analysisDir, sampleMetadata$sampleID, '.Tn5_insertion_points.tagAlign.gz', sep="")
sampleMetadata        <- mutate(sampleMetadata, tn5insertionPointsBed=tn5insertionPointsBed)
filtMetadataOrdered   <- sampleMetadata
conditionsToCount     <- c("EtOH-nlDensity","RA-high", "TGFb-high")
nConditions           <- length(conditionsToCount)
nReplicates           <- 3

countMatrix           <- as.matrix(counts(fragmentCounts))
colnames(countMatrix) <- str_extract(colnames(fragmentCounts), pattern = "[^\\.]+")
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
