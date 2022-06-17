tss <- read_tsv(here('refs', 'EnsgToHg38CanonicalTssMapping.tsv'))
bigTib <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))

tss.bed <- tss %>% 
  dplyr::select(chrom, tx_start, tx_end)

write_tsv(tss.bed, file = here('refs', 'tss.bed'), col_names = FALSE)

allPeaks <- bigTib %>%
  dplyr::select(peak) %>% 
  tidyr::separate(col = peak, into = c("chrom", "start", "end")) %>% 
  drop_na()

write_tsv(allPeaks, here('extractedData','allPeaks.bed'), col_names = FALSE)
