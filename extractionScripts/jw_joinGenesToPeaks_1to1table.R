cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksAnno     <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'jw_peakAnno_AtacPeaksMergedSummitWindows_mergedist50.tsv'))
  peaksVal      <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'ATACseq', 'jw_AtacPeaksMergedSummitWindows_mergedist50.tsv'))
  deseqTib      <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'RNAseq', 'IL1a_1hr_KB_normalizedCounts_meltedData.tsv'))
  outputLoc     <- here('extractedData', 'Jurida_Weiterer', 'jw_genePeak_1to1mapping_mergedDist50.tsv')
} else {
  peaksAnno     <- cmdargs[1]
  peaksTib      <- cmdargs[2]
  deseqTib      <- cmdargs[3]
  outputLoc     <- cmdargs[4]
}

#combine peak annotations and peak values
peaksAnno <- peaksAnno %>% 
  dplyr::mutate(peak = paste(seqnames, start, end, sep=":"))

peaksTib <- dplyr::left_join(peaksVal, peaksAnno)

peaksTib <- peaksTib %>% 
  dplyr::rename(ensg = ENSEMBL)

deseqTib <- deseqTib %>% 
  dplyr::select(gene_id, sampleID, tpm) %>% 
  tidyr::pivot_wider(names_from = sampleID, values_from = tpm) %>% 
  dplyr::rename(untreated_rna = control, IL1a_1hr_rna = IL1a_1hr, ENSEMBL = gene_id)

fullTib <- dplyr::left_join(deseqTib, peaksTib) %>% 
  dplyr::rename(untreated_normATAC = untreated, IL1a_1hr_normATAC = il1a_1hr)

write_tsv(fullTib, file = outputLoc)
