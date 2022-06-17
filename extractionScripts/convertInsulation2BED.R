cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  tadFile <- read_tsv(here('extractedData', 'Hi-C_MCF7_MCF10A_processed_HiCfiles', 'TAD_boundaries', 'HiCStein-MCF7_allChr_hg19.insulation.boundaries'), col_names = FALSE)
  outputFile <- here('extractedData', 'HiCStein-MCF7_allChr_hg19.bed')
} else {
  tadFile <- read_tsv(cmdargs[1])  
  outputFile <- cmdargs[2]
}

colnames(tadFile) <- c('header', 'start', 'end', 'binStart', 'binEnd', 'binMidpoint', 'header2', 'insulationScore')

tadFile <- tadFile %>% 
  tidyr::separate(col = header, sep = "\\|", into = c(NA, NA, 'tad'))

tadFile <- tadFile %>% 
  tidyr::separate(tad, into = c("chr", "start", "end")) %>% 
  dplyr::select(chr, start, end)


write_tsv(tadFile, file = outputFile, col_names = FALSE)
