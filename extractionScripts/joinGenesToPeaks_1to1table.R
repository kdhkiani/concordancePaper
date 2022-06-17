cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksAnno     <- read_tsv(here('extractedData', 'peakAnno_allCondsMergedSummitWindows_mergedist50.tsv'))
  peaksVal      <- read_tsv(here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'))
  deseqTib      <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
  outputLoc     <- here('extractedData', 'genePeak_1to1mapping_mergedDist50.tsv')
} else {
  peaksAnno     <- cmdargs[1]
  peaksTib      <- cmdargs[2]
  deseqTib      <- cmdargs[3]
  outputLoc     <- cmdargs[4]
}

#combine peak annotations and peak values
peaksAnno <- peaksAnno %>% 
  dplyr::mutate(peak = paste(seqnames, start, end, sep=":"))

peaksVal <- peaksVal %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":"))

peaksTib <- dplyr::left_join(peaksVal, peaksAnno)

peaksTib <- peaksTib %>% 
  dplyr::rename(ensg = ENSEMBL)

fullTib <- dplyr::left_join(deseqTib, peaksTib, by = "ensg")

write_tsv(fullTib, file = outputLoc)
