cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peakTib      <- read_tsv(here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.annotated.tsv'))
  diffpeaks    <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.tsv'))
  outputFolder <- here('extractedData', 'HINT', '')
 
} else {
  fragmentCounts <- readRDS(cmdargs[1])  
  outputFile <- cmdargs[2]
}

diffpeaks <- diffpeaks %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":"))

RA.peaks <- diffpeaks %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE) %>% 
  dplyr::select(peak)

RA.peaks <- RA.peaks$peak

TGFb.peaks <- diffpeaks %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE) %>% 
  dplyr::select(peak)

TGFb.peaks <- TGFb.peaks$peak

peakTib <- peakTib %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":")) %>% 
  dplyr::mutate(RA.filt = peak %in% RA.peaks) %>% 
  dplyr::mutate(TGFb.filt = peak %in% TGFb.peaks)

RA.out <- peakTib %>% 
  dplyr::filter(RA.filt == FALSE) %>% 
  dplyr::filter(`RA-high-avgNormFragmentCounts` > 30) %>% 
  dplyr::select(chrom, startLocs, endLocs)

write_tsv(RA.out, file = paste0(outputFolder, 'RAhigh.bed'), col_names = FALSE)

TGFb.out <- peakTib %>% 
  dplyr::filter(TGFb.filt == FALSE) %>% 
  dplyr::filter(`TGFb-high-avgNormFragmentCounts` > 30) %>% 
  dplyr::select(chrom, startLocs, endLocs)

write_tsv(TGFb.out, file = paste0(outputFolder, 'TGFbhigh.bed'), col_names = FALSE)

EtOH.out <- peakTib %>% 
  dplyr::filter(`EtOH-nlDensity-avgNormFragmentCounts` > 30) %>% 
  dplyr::select(chrom, startLocs, endLocs)

write_tsv(EtOH.out, file = paste0(outputFolder, 'EtOH.bed'), col_names = FALSE)


