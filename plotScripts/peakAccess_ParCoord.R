theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigTib        <- fread(here('extractedData', '50KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                         select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc",
                                    "EtOH-nlDensity-avgNormFragmentCounts", "RA-high-avgNormFragmentCounts",
                                    "RA-med-avgNormFragmentCounts", "RA-low-avgNormFragmentCounts",
                                    "TGFb-high-avgNormFragmentCounts", "TGFb-med-avgNormFragmentCounts",
                                    "TGFb-low-avgNormFragmentCounts")) #read in only specific columns to avoid choking up memory
  diffPeaks     <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.tsv'))
  outputFolder   <- here("plots", "validation")
} else {
  geneTib        <- read_tsv(cmdargs[1])
  samplemetadata <- read_tsv(cmdargs[2])
}

#convert to tibble
bigTib <- as_tibble(bigTib)

#Add unique peak identifier
bigTib <- mutate(bigTib, peak_loc = paste(peak_chrom, peak_startLoc, peak_endLoc,sep = ":"))

bigTib <- bigTib %>% 
    rowwise() %>% #to do minimum across rows instead of columns
  dplyr::mutate(RA.maxNormFragmentCounts = max(`EtOH-nlDensity-avgNormFragmentCounts`,`RA-high-avgNormFragmentCounts`,
                                            `RA-med-avgNormFragmentCounts`, `RA-low-avgNormFragmentCounts`)) 
##Plot Par Coord for RA
gene <- "SLC5A5"

plotTib <- bigTib %>%
  filter(gene_name == gene) %>% 
  dplyr::select(matches("(gene_name)|(avgNormFragmentCounts)|(peak_loc)|(maxNormFragmentCounts)")) %>% 
  dplyr::select(-contains("TGFb"))

##filter
plotTib <- plotTib %>% dplyr::filter(RA.maxNormFragmentCounts > 30)


#change column order
plotTib <- plotTib %>% dplyr::relocate(`RA-high-avgNormFragmentCounts`, .after = `RA-med-avgNormFragmentCounts`) %>% 
  dplyr::relocate(`RA-low-avgNormFragmentCounts`, .before = `RA-med-avgNormFragmentCounts`)

#annotate differential peaks
RA.diffPeaks <- diffPeaks %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE) %>% 
  dplyr::mutate(peak_loc = paste(chrom, startLocs, endLocs, sep=":"))

plotTib <- plotTib %>% 
  dplyr::mutate(isDiffPeak = peak_loc %in% RA.diffPeaks$peak_loc)

#plot raw values
p1 <- ggparcoord(plotTib, columns = 2:5, groupColumn = 8, alphaLines = 1, scale = "globalminmax", showPoints = T) +
  theme(legend.position = "none") +
  xlab("") + ylab("Accessibility (normalized fragment counts)")

outputLoc <- paste0(outputFolder,"/","SLC5A5_peakAccessibility.svg")
ggsave(outputLoc, p1)

##Plot Par Coord for second gene
gene <- "HOXA1"

plotTib <- bigTib %>%
  filter(gene_name == gene) %>% 
  dplyr::select(matches("(gene_name)|(avgNormFragmentCounts)|(peak_loc)|(maxNormFragmentCounts)")) %>% 
  dplyr::select(-contains("TGFb"))

##filter
plotTib <- plotTib %>% dplyr::filter(RA.maxNormFragmentCounts > 30)
#change column order
plotTib <- plotTib %>% dplyr::relocate(`RA-low-avgNormFragmentCounts`, .after = `EtOH-nlDensity-avgNormFragmentCounts`) %>% 
  dplyr::relocate(`RA-high-avgNormFragmentCounts`, .after = `RA-med-avgNormFragmentCounts`)

#annotate differential peaks
RA.diffPeaks <- diffPeaks %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE) %>% 
  dplyr::mutate(peak_loc = paste(chrom, startLocs, endLocs, sep=":"))

plotTib <- plotTib %>% 
  dplyr::mutate(isDiffPeak = peak_loc %in% RA.diffPeaks$peak_loc)

#plot
p2 <- ggparcoord(plotTib, columns = 2:5, groupColumn = 8, alphaLines = 1, scale = "globalminmax", showPoints = T) +
  theme(legend.position = "none") +
  xlab("") + ylab("Accessibility (normalized fragment counts)")

outputLoc <- paste0(outputFolder,"/","HOXA1_peakAccessibility.svg")
ggsave(outputLoc, p2)
