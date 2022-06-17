theme_set(theme_classic())

########
#using source data from paper website
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib       <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'jw_genePeak_1to1mapping_mergedDist50.tsv'))
  outputFolder <- here('plots', 'proportion')
} else {
  bigtib         <- cmdargs[1]
  outputFolder   <- cmdargs[2]
}

##massage and prepare data from bigtib

#drop any peaks with no genes assigned
bigtib <- drop_na(bigtib)

#add 1 to avoid infinite values when doing log2 fold changes
bigtib <- bigtib %>% 
  dplyr::mutate(untreated_rna = untreated_rna + 1) %>% 
  dplyr::mutate(IL1a_1hr_rna = IL1a_1hr_rna + 1) %>% 
  dplyr::mutate(untreated_normATAC = untreated_normATAC + 1) %>% 
  dplyr::mutate(IL1a_1hr_normATAC = IL1a_1hr_normATAC + 1)

#annotate max peak value between control and IL1alpha treated conditions
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(maxValATAC = max(untreated_normATAC, IL1a_1hr_normATAC)) %>% 
  dplyr::mutate(maxValRNA = max(untreated_rna, IL1a_1hr_rna))

#create histogram of distribution of max value between two samples of RNA
p1 <- ggplot(bigtib, aes(x = maxValATAC)) + 
  geom_histogram(fill = "light gray", color = "black", alpha = 0.8) + 
  xlim(-1,50) + geom_vline(xintercept = 5)

#outputLoc <- paste0(outputFolder, "/maxRNAhist.svg")
#ggsave(plot = p1, filename = outputLoc, height = 5, width = 5, units = "in", dpi = 300)

#remove genes where max rna between cell populations is less than 5 tpm
bigtib <- bigtib %>% 
  dplyr::filter(maxValRNA > 5)

#calculate log2FC for accessibility and create a filtering column
bigtib <- bigtib %>% 
  dplyr::mutate(log2FCatac = log2(IL1a_1hr_normATAC/untreated_normATAC)) %>% 
  dplyr::mutate(peakCutoff = abs(log2FCatac) > 0.58)

#annotate peak change direction
bigtib <- bigtib %>% 
  dplyr::mutate(peakDir = case_when(IL1a_1hr_normATAC > untreated_normATAC ~ "UP",
                                    TRUE ~ "DOWN"))
#determine log2FC for expression
bigtib <- bigtib %>% 
  dplyr::mutate(log2FC = log2(IL1a_1hr_rna/untreated_rna))

#arrange based on log2FC expression values
bigtib <- bigtib %>% 
  dplyr::arrange(desc(log2FC))

#count number of peak combinations
bigtib <- bigtib %>% 
  group_by(ENSEMBL) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & peakCutoff == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & peakCutoff == FALSE))

#get counts for number of peaks per gene
peakCounts <- bigtib %>% dplyr::count(ENSEMBL)

#get counts of all peaks per gene and sig peaks
plotTib <- bigtib %>% 
  dplyr::select(ENSEMBL, log2FC, Up_Sig, Up_NotSig, Down_Sig, Down_NotSig) %>% unique()

#change to proportions
plotTib <- plotTib %>% 
  dplyr::mutate(totalPeaks = sum(Up_Sig, Up_NotSig, Down_Sig, Down_NotSig)) %>% 
  
  dplyr::mutate(Up_Sig = Up_Sig/totalPeaks) %>% 
  dplyr::mutate(Up_NotSig = Up_NotSig/totalPeaks) %>% 
  dplyr::mutate(Down_Sig = Down_Sig/totalPeaks) %>% 
  dplyr::mutate(Down_NotSig = Down_NotSig/totalPeaks) %>% 
  
  dplyr::select(- totalPeaks)

#annotate gene complexity type based on defintions from gonzalez et al, 2015
plotTib <- left_join(plotTib, peakCounts) %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::mutate(lociComplexity = case_when(numPeaks > 7 ~ "high",
                                           TRUE ~ "low"))

#change table to to long format

plotTib <- plotTib %>% gather("peakType", "proportion", -ENSEMBL, -log2FC, -numPeaks, -lociComplexity)


#change plotTib to attempt facetting 
plotTib <- plotTib %>% 
  tidyr::separate(peakType, into = c("direction", "significantPeak"))


#for clarity, filter out non significant peaks

Sig.plotTib <- plotTib %>% 
  dplyr::filter(significantPeak == "Sig") 

### get rank values
top50 <- Sig.plotTib %>%
  arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(1:50) %>%
  dplyr::select(ENSEMBL, log2FC, lociComplexity) %>% unique()
top50 <- top50$ENSEMBL

top100 <- Sig.plotTib %>%
  arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(51:100) %>%
  dplyr::select(ENSEMBL, log2FC, lociComplexity) %>% unique()
top100 <- top100$ENSEMBL

bot50 <- Sig.plotTib %>%
  arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(1:50) %>%
  dplyr::select(ENSEMBL, log2FC, lociComplexity) %>% unique()
bot50 <- bot50$ENSEMBL

bot100<- Sig.plotTib %>%
  arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(51:100) %>%
  dplyr::select(ENSEMBL, log2FC, lociComplexity) %>% unique()
bot100 <- bot100$ENSEMBL

##annotate based on rank position
Sig.plotTib <- Sig.plotTib %>% 
  dplyr::mutate(geneRank = case_when(ENSEMBL %in% top50 ~ "Top 50 IL1a stim",
                                     ENSEMBL %in% top100 ~ "Next 50 IL1a stim",
                                     ENSEMBL %in% bot50 ~ "Top 50 untreated",
                                     ENSEMBL %in% bot100 ~ "Next 50 untreated",
                                     TRUE ~ "other"))

Sig.plotTib <- Sig.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 untreated", "Top 50 untreated",
                                                       "other",
                                                       "Top 50 IL1a stim", "Next 50 IL1a stim"
                                                       )))
hi.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter(lociComplexity == "high")

p2 <- ggscatter(hi.Sig.plotTib, x = "log2FC", y = "proportion",
                facet.by = c("direction"), color = "geneRank",
                palette = c("#CACECB", "#808782", "#000000", "#EFC22E", "#F8E4A0"),
                alpha = 0.8, size = 1.5) +
  ggtitle("Proportion of ATAC peaks differenitally accessible for top 100 untreated and IL1alpha treated genes from Weiterer et al., 2020") +
  xlab("log2 Fold Change in Expression") +
  ylab("Proportion of DHS peaks differentially accessible") +
  ylim(0,1)

ggsave(p2, filename = paste0(outputFolder,"/jwProportionPlot_highComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)

lo.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter(lociComplexity == "low")

p3 <- ggscatter(lo.Sig.plotTib, x = "log2FC", y = "proportion",
                facet.by = c("direction"), color = "geneRank",
                palette = c("#CACECB", "#808782", "#000000", "#EFC22E", "#F8E4A0"),
                alpha = 0.8, size = 1.5) +
  ggtitle("Proportion of DHS peaks differenitally accessible for top 100 CD34 and CD14 genes from Gonzalez et al., 2015") +
  xlab("log2 Fold Change in Expression") +
  ylab("Proportion of DHS peaks differentially accessible") +
  ylim(0,1)

ggsave(p3, filename = paste0(outputFolder,"/jwProportionPlot_lowComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Up" & lociComplexity == "high") %>% 
  dplyr::filter(geneRank == "Top 50 IL1a stim" | geneRank == "Next 50 IL1a stim")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Down" & lociComplexity == "high") %>% 
  dplyr::filter(geneRank == "Top 50 untreated" | geneRank == "Next 50 untreated")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank) %>% 
  bind_rows(hist.plotTib)

p4 <- hist.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 untreated", "Top 50 untreated","Next 50 IL1a stim", "Top 50 IL1a stim"))) %>% 
  ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
  geom_density_ridges(alpha = 0.8, quantile_lines = T, quantiles=2) +
  scale_fill_manual(values = c(c("#CACECB", "#808782", "#F8E4A0", "#EFC22E"))) +
  theme(legend.position = "")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))

ggsave(p4, filename = paste0(outputFolder,"/jwProportionDensity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)


tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Up" & lociComplexity == "low") %>% 
  dplyr::filter(geneRank == "Top 50 IL1a stim" | geneRank == "Next 50 IL1a stim")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Down" & lociComplexity == "low") %>% 
  dplyr::filter(geneRank == "Top 50 untreated" | geneRank == "Next 50 untreated")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank) %>% 
  bind_rows(hist.plotTib)

p5 <- hist.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 untreated", "Top 50 untreated","Next 50 IL1a stim", "Top 50 IL1a stim"))) %>% 
  ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
  geom_density_ridges(alpha = 0.8, quantile_lines = T, quantiles=2) +
  scale_fill_manual(values = c(c("#CACECB", "#808782", "#F8E4A0", "#EFC22E"))) +
  theme(legend.position = "")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))

ggsave(p5, filename = paste0(outputFolder,"/jwProportionDensity_lowComplexity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
