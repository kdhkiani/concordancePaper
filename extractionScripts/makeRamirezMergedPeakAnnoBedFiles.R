##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peakstib   <- read_tsv(here('extractedData', 'Ramirez', 'ATACseq', 'ramirez_AtacPeaksMergedSummitWindows_mergedist50.tsv'))
  outputLoc <- here('extractedData', 'Ramirez', 'ATACseq')
} else {
  peaksanno    <- read_tsv(cmdargs[1])
  outputLoc  <- cmdargs[2]
}

#pull out the relevant values 
peakstib <- peakstib %>% 
  dplyr::select(matches("(peak)|(HL60)|(monomac)")) 

#annotate max value between RA/etOH and TGFb/etOH
peakstib <- peakstib %>% 
  rowwise() %>% 
  dplyr::mutate(maxVal = max(across(where(is.double)))) 
  
#filter

peakstib <- peakstib %>% 
  dplyr::filter(maxVal > 5)

###create bed files for peaks
bed_file_fields <- peakstib %>% 
  tidyr::separate(peak, into = c("chrom", "startLocs", "endLocs")) %>% 
  dplyr::select(chrom, startLocs, endLocs) 
write_tsv(bed_file_fields, paste0(outputLoc,"/HL60_monomac_peaks.bed"), col_names = F)
