window40_tib <- read_tsv(here('extractedData','40KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

window100_tib <- read_tsv(here('extractedData','100KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

out <- window100_tib %>% anti_join(window40_tib)

write_tsv(x = out, file = here('extractedData', 'distalWindow_mergedist50_joinedTablePeaksNearGenes.tsv'))
