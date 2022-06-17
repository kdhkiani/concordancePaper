##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peakstib   <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.bed'))
  outputLoc <- here('extractedData')
} else {
  peaksanno    <- read_tsv(cmdargs[1])
  outputLoc  <- cmdargs[2]
}

#pull out the relevant values 
peakstib <- peakstib %>% 
  dplyr::select(matches('(chrom)|(startLocs)|(endLocs)|(RA-high-isDiffPeak)|(TGFb-high-isDiffPeak)|(high-avgNormFragmentCounts)|(nlDensity-avgNormFragmentCounts)')) %>% 
  dplyr::select(- matches('and'))

#annotate max value between RA/etOH and TGFb/etOH
peakstib <- peakstib %>% 
  rowwise() %>% 
  dplyr::mutate(RA.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`RA-high-avgNormFragmentCounts`)) %>% 
  dplyr::mutate(TGFb.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`TGFb-high-avgNormFragmentCounts`))

###create bed files for diff peaks

##no filter

#RA high peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/RAhigh_diffPeak_maxVal30.bed"), col_names = F)

#TGFb peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/TGFbhigh_diffPeak_maxVal30.bed"), col_names = F)

##maxVal > 50

#RA high peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE) %>% 
  dplyr::filter(RA.maxVal > 50) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/RAhigh_diffPeak_maxVal50.bed"), col_names = F)

#TGFb peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE) %>% 
  dplyr::filter(TGFb.maxVal > 50) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/TGFbhigh_diffPeak_maxVal50.bed"), col_names = F)
