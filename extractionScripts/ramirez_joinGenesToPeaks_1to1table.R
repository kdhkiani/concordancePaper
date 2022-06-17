cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksAnno     <- read_tsv(here('extractedData', 'Ramirez', 'ATACseq', 'ramirez_peakAnno_AtacPeaksMergedSummitWindows_mergedist50.tsv'))
  peaksVal      <- read_tsv(here('extractedData', 'Ramirez', 'ATACseq', 'ramirez_AtacPeaksMergedSummitWindows_mergedist50.tsv'))
  deseqTib      <- read_tsv(here('extractedData', 'Ramirez', 'RNAseq', 'HL60_normCounts_meltedData.tsv'))
  outputLoc     <- here('extractedData', 'Ramirez', 'ramirez_genePeak_1to1mapping_mergedDist50.tsv')
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

#add some more specificity to column names
colnames(peaksTib)[2:16] <- paste0(colnames(peaksTib)[2:16],"_normFragmentCounts")

deseqTib <- deseqTib %>% 
  dplyr::select(gene_id, sampleID, tpm) %>% 
  tidyr::pivot_wider(names_from = sampleID, values_from = tpm) %>% 
  dplyr::rename(ensg = gene_id)

##rename columns so they make more sense
deseqTib <- deseqTib %>% 
  dplyr::rename(SRR3214253_HL60_tpm=SRR3214253) %>%
  dplyr::rename(SRR3214254_HL60_tpm=SRR3214254) %>%
  dplyr::rename(SRR3214255_HL60_tpm=SRR3214255) %>%
  
  dplyr::rename(SRR3214325_monocyte_tpm=SRR3214325) %>%
  dplyr::rename(SRR3214326_monocyte_tpm=SRR3214326) %>%
  dplyr::rename(SRR3214327_monocyte_tpm=SRR3214327) %>%

  dplyr::rename(SRR3214289_neutrophil12hr_tpm=SRR3214289) %>%
  dplyr::rename(SRR3214290_neutrophil12hr_tpm=SRR3214290) %>%
  dplyr::rename(SRR3214291_neutrophil12hr_tpm=SRR3214291) %>%
  
  dplyr::rename(SRR3214301_neutrophil120hr_tpm=SRR3214301) %>%
  dplyr::rename(SRR3214302_neutrophil120hr_tpm=SRR3214302) %>%
  dplyr::rename(SRR3214303_neutrophil120hr_tpm=SRR3214303) %>%
  
  dplyr::rename(SRR3214343_monomac_tpm=SRR3214343) %>%
  dplyr::rename(SRR3214344_monomac_tpm=SRR3214344) %>%
  dplyr::rename(SRR3214345_monomac_tpm=SRR3214345)

fullTib <- dplyr::left_join(deseqTib, peaksTib)

write_tsv(fullTib, file = outputLoc)
