##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peakstib   <- read_tsv(here('extractedData', 'allCondsMergedSummitWindows_mergedist50.bed'))
  outputLoc <- here('extractedData')
} else {
  peaksanno    <- read_tsv(cmdargs[1])
  outputLoc  <- cmdargs[2]
}

#pull out the relevant values 
peakstib <- peakstib %>% 
  dplyr::select(matches('(chrom)|(startLocs)|(endLocs)|(high-avgNormFragmentCounts)|(nlDensity-avgNormFragmentCounts)')) %>% 
  dplyr::select(- matches('and'))

#annotate max value between RA/etOH and TGFb/etOH
peakstib <- peakstib %>% 
  rowwise() %>% 
  dplyr::mutate(RA.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`RA-high-avgNormFragmentCounts`)) %>% 
  dplyr::mutate(TGFb.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`TGFb-high-avgNormFragmentCounts`))

###create bed files for peaks
bed_file_fields <- peakstib %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/allPeaks_mergedPeaks_noFilter.bed"), col_names = F)


##maxVal > 30

#RA high peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(RA.maxVal > 30) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/RAhigh_mergedPeaks_maxVal30.bed"), col_names = F)

#TGFb peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(TGFb.maxVal > 30) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/TGFbhigh_mergedPeaks_maxVal30.bed"), col_names = F)

#EtOH normal density peaks
bed_file_fields <- peakstib %>% 
  dplyr::filter(`EtOH-nlDensity-avgNormFragmentCounts` > 30) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/EtOHnlDensity_mergedPeaks_maxVal30.bed"), col_names = F)
