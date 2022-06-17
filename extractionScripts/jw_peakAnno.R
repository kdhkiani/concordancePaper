cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksFile     <- here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'jw_AtacPeaksMergedSummitWindows_mergedist50.bed')
  outputLoc     <- here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'jw_peakAnno_AtacPeaksMergedSummitWindows_mergedist50.tsv')
} else {
  peaksFile     <- cmdargs[1]
  outputLoc     <- cmdargs[2]
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnno <- annotatePeak(peaksFile, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnnoTib <- as_tibble(peakAnno@anno)

#adjust start and end locs back 1 bp so sites will match up given 0/1 based indexing discrepancies
peakAnnoTib <- peakAnnoTib %>% 
  dplyr::mutate(start = start -1) %>% 
  dplyr::mutate(end = end)

write_tsv(x = peakAnnoTib, file = outputLoc)
