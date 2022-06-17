cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksAnno     <- read_tsv(here('extractedData', 'peakAnno_allCondsMergedSummitWindows_mergedist50.tsv'))
  peaksVal      <- read_tsv(here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.annotated.tsv'))
  deseqTib      <- read_tsv(here('extractedData', 'atac_intron_meltedData.tsv'))
  outputLoc     <- here('extractedData', 'intron_genePeak_1to1mapping_mergedDist50.tsv')
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

#make intron deseqtib wide
deseqTib <- deseqTib %>% 
  dplyr::select(gene_id, sampleID, counts) %>% 
  tidyr::pivot_wider(names_from = sampleID, values_from = counts) %>% 
  rename(ensg = gene_id)

colnames(deseqTib)[-1] = paste0(colnames(deseqTib)[-1],"_intron")

fullTib <- dplyr::left_join(deseqTib, peaksTib)

write_tsv(fullTib, file = outputLoc)
