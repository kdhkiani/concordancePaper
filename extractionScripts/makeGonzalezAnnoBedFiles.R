##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksanno   <- read_tsv(here('Gonzalez2015', 'peaksTable.tsv'))
  peaksdat    <- read_tsv(here('Gonzalez2015', 'DNaseCnts.tsv'))
  outputLoc <- here('extractedData')
} else {
  peaksanno    <- read_tsv(cmdargs[1])
  peaksdat <- read_tsv(cmdargs[2])
  outputLoc  <- cmdargs[3]
}

#add annotation for peak types in peaksanno table
peaksanno <- peaksanno %>% 
  dplyr::mutate(isCD34peak = grepl(x=accessPattern, pattern = "CD34")) %>% 
  dplyr::mutate(isCD14peak = grepl(x=accessPattern, pattern = "CD14"))

#calculate log2FC for peaks and whether peak passes various cutoffs
peaksdat <- peaksdat %>% 
  dplyr::mutate(log2fc = log2(CD14/CD34)) %>% 
  dplyr::mutate(fcCutoff = abs(log2fc) > 1.5) %>% 
  rowwise() %>% 
  dplyr::mutate(maxVal = max(CD34,CD14))

#join tables together into master table
peaktib <- left_join(peaksanno, peaksdat, by = "peakID")
                
###create bed files for peak calls

##CD34 peaks
bed_file_fields <- peaktib %>% 
  dplyr::filter(isCD34peak == TRUE) %>% 
  dplyr::select(chr, start, end) 
write_tsv(bed_file_fields, paste0(outputLoc,"/CD34.bed"), col_names = F)

##CD14 peaks
bed_file_fields <- peaktib %>% 
  dplyr::filter(isCD14peak == TRUE) %>% 
  dplyr::select(chr, start, end) 
write_tsv(bed_file_fields, paste0(outputLoc,"/CD14.bed"), col_names = F)

##diffpeaks based on |log2fc| > 2, no filter
bed_file_fields <- peaktib %>% 
  dplyr::filter(fcCutoff == TRUE) %>% 
  dplyr::select(chr, start, end) 
write_tsv(bed_file_fields, paste0(outputLoc,"/gonzalezDiffPeak_noFilter.bed"), col_names = F)

#diffpeaks based on |log2fc| > 2, max value > 10
bed_file_fields <- peaktib %>% 
  dplyr::filter(fcCutoff == TRUE) %>% 
  dplyr::filter(maxVal > 10) %>% 
  dplyr::select(chr, start, end) 
write_tsv(bed_file_fields, paste0(outputLoc,"/gonzalezDiffPeak_10Filter.bed"), col_names = F)

#diffpeaks based on |log2fc| > 2, max value > 30
bed_file_fields <- peaktib %>% 
  dplyr::filter(fcCutoff == TRUE) %>% 
  dplyr::filter(maxVal > 30) %>% 
  dplyr::select(chr, start, end) 
write_tsv(bed_file_fields, paste0(outputLoc,"/gonzalezDiffPeak_30Filter.bed"), col_names = F)

#diffpeaks based on |log2fc| > 2, max value > 50
bed_file_fields <- peaktib %>% 
  dplyr::filter(fcCutoff == TRUE) %>% 
  dplyr::filter(maxVal > 50) %>% 
  dplyr::select(chr, start, end) 
write_tsv(bed_file_fields, paste0(outputLoc,"/gonzalezDiffPeak_50Filter.bed"), col_names = F)
