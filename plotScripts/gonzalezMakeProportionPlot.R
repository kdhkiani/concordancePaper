theme_set(theme_classic())

########
#using source data from paper website
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksanno   <- read_tsv(here('extractedData', 'Gonzalez2015', 'peaksTable.tsv'))
  peaksdat    <- read_tsv(here('extractedData', 'Gonzalez2015', 'DNaseCnts.tsv'))
  deseqtib    <- read_tsv(here('extractedData', 'Gonzalez2015', 'RNAseqCnts.tsv')) 
  outputFolder <- here('plots', 'proportion')
} else {
  peaksanno   <- read_tsv(cmdargs[1])
  peaksdat    <- read_tsv(cmdargs[2])
  deseqtib    <- read_tsv(cmdargs[3])
  outputFolder   <- cmdargs[4]
}

peaktib <- left_join(peaksanno, peaksdat)

peaktib <- peaktib %>% 
  dplyr::rename(CD34_dhs = CD34, CD14_dhs= CD14) %>% 
  dplyr::select(- H1hesc, - CD56, - CD3, -CD19)

deseqtib <- deseqtib %>% 
  dplyr::rename(CD34_rna = CD34, CD14_rna= CD14) %>% 
  dplyr::select(- H1hesc, - CD56, - CD3, -CD19)

bigtib <- left_join(peaktib, deseqtib, by = "symbol")

bigtib <- bigtib %>% 
  dplyr::mutate(peakFilter = grepl(pattern = "CD[1-3]4", x = accessPattern)) %>% 
  dplyr::filter(peakFilter == TRUE)

##massage and prepare data from bigtib

#add 1 to avoid infinite values when doing log2 fold changes
bigtib <- bigtib %>% 
  dplyr::mutate(CD34_rna = CD34_rna + 1) %>% 
  dplyr::mutate(CD34_dhs = CD34_dhs + 1) %>% 
  dplyr::mutate(CD14_rna = CD14_rna + 1) %>% 
  dplyr::mutate(CD14_dhs = CD14_dhs + 1)

#annotate max peak value between CD34, CD14 
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(maxValDHS = max(CD14_dhs, CD34_dhs)) %>% 
  dplyr::mutate(maxValRNA = max(CD14_rna, CD34_rna))

#create histogram of distribution of max value between two samples of RNA
p1 <- ggplot(bigtib, aes(x = maxValRNA)) + 
  geom_histogram(fill = "light gray", color = "black", alpha = 0.8) + 
  xlim(-1,50) + geom_vline(xintercept = 5)

outputLoc <- paste0(outputFolder, "/maxRNAhist.svg")
ggsave(plot = p1, filename = outputLoc, height = 5, width = 5, units = "in", dpi = 300)

#remove genes where max rna between cell populations is less than 5 tpm
bigtib <- bigtib %>% 
  dplyr::filter(maxValRNA > 5)

#calculate log2FC for accessibility and create a filtering column
bigtib <- bigtib %>% 
  dplyr::mutate(log2FCdhs = log2(CD14_dhs/CD34_dhs)) %>% 
  dplyr::mutate(peakCutoff = abs(log2FCdhs) > 2)

#annotate peak change direction
bigtib <- bigtib %>% 
  dplyr::mutate(peakDir = case_when(CD14_dhs > CD34_dhs ~ "UP",
                                    TRUE ~ "DOWN"))
#determine log2FC for expression
bigtib <- bigtib %>% 
  dplyr::mutate(log2FC = log2(CD14_rna/CD34_rna))

#arrange based on log2FC expression values
bigtib <- bigtib %>% 
  dplyr::arrange(desc(log2FC))

#count number of peak combinations
bigtib <- bigtib %>% 
  group_by(symbol) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & peakCutoff == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & peakCutoff == FALSE))

#get counts for number of peaks per gene
peakCounts <- bigtib %>% dplyr::count(symbol)

#get counts of all peaks per gene and sig peaks
plotTib <- bigtib %>% 
  dplyr::select(symbol, log2FC, Up_Sig, Up_NotSig, Down_Sig, Down_NotSig) %>% unique()

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

plotTib <- plotTib %>% gather("peakType", "proportion", -symbol, -log2FC, -numPeaks, -lociComplexity)


#change plotTib to attempt facetting 
plotTib <- plotTib %>% 
  tidyr::separate(peakType, into = c("direction", "significantPeak"))


#for clarity, filter out non significant peaks

Sig.plotTib <- plotTib %>% 
  dplyr::filter(significantPeak == "Sig") 

### get rank values
top50 <- Sig.plotTib %>%
  arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>%
  dplyr::select(symbol, log2FC, lociComplexity) %>% unique()
top50 <- top50$symbol

top100 <- Sig.plotTib %>%
  arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>%
  dplyr::select(symbol, log2FC, lociComplexity) %>% unique()
top100 <- top100$symbol

bot50 <- Sig.plotTib %>%
  arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>%
  dplyr::select(symbol, log2FC, lociComplexity) %>% unique()
bot50 <- bot50$symbol

bot100<- Sig.plotTib %>%
  arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>%
  dplyr::select(symbol, log2FC, lociComplexity) %>% unique()
bot100 <- bot100$symbol

##annotate based on rank position
Sig.plotTib <- Sig.plotTib %>% 
  dplyr::mutate(geneRank = case_when(symbol %in% top50 ~ "Top 50 CD14",
                                     symbol %in% top100 ~ "Next 50 CD14",
                                     symbol %in% bot50 ~ "Top 50 CD34",
                                     symbol %in% bot100 ~ "Next 50 CD34",
                                     TRUE ~ "other"))

Sig.plotTib <- Sig.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Top 50 CD34", "Next 50 CD34",
                                                       "other",
                                                       "Next 50 CD14", "Top 50 CD14")))
hi.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter(lociComplexity == "high")

p2 <- ggscatter(hi.Sig.plotTib, x = "log2FC", y = "proportion",
                facet.by = c("direction"), color = "geneRank",
                palette = c("#3358FF", "#ADBEFF", "#000000","#F9D3B4", "#EE8833"),
                alpha = 0.8, size = 1.5) +
  ggtitle("Proportion of DHS peaks differenitally accessible for top 100 CD34 and CD14 genes from Gonzalez et al., 2015") +
  xlab("log2 Fold Change in Expression") +
  ylab("Proportion of DHS peaks differentially accessible")

ggsave(p2, filename = paste0(outputFolder,"/gonzalezProportionPlot_highComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)

lo.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter(lociComplexity == "low")

p3 <- ggscatter(lo.Sig.plotTib, x = "log2FC", y = "proportion",
                facet.by = c("direction"), color = "geneRank",
                palette = c("#3358FF", "#ADBEFF", "#000000","#F9D3B4", "#EE8833"),
                alpha = 0.8, size = 1.5) +
  ggtitle("Proportion of DHS peaks differenitally accessible for top 100 CD34 and CD14 genes from Gonzalez et al., 2015") +
  xlab("log2 Fold Change in Expression") +
  ylab("Proportion of DHS peaks differentially accessible")

ggsave(p3, filename = paste0(outputFolder,"/gonzalezProportionPlot_lowComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Up" & lociComplexity == "high") %>% 
  dplyr::filter(geneRank == "Top 50 CD14" | geneRank == "Next 50 CD14")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Down" & lociComplexity == "high") %>% 
  dplyr::filter(geneRank == "Top 50 CD34" | geneRank == "Next 50 CD34")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank) %>% 
  bind_rows(hist.plotTib)

p4 <- hist.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 CD34", "Top 50 CD34", "Next 50 CD14", "Top 50 CD14"))) %>% 
  ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
  geom_density_ridges(alpha = 0.8, quantile_lines = T, quantiles=2) +
  scale_fill_manual(values = c(c("#ADBEFF", "#3358FF", "#F9D3B4", "#EE8833"))) +
  theme(legend.position = "")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))

ggsave(p4, filename = paste0(outputFolder,"/gonzalezProportionDensity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)


tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Up" & lociComplexity == "low") %>% 
  dplyr::filter(geneRank == "Top 50 CD14" | geneRank == "Next 50 CD14")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Down" & lociComplexity == "low") %>% 
  dplyr::filter(geneRank == "Top 50 CD34" | geneRank == "Next 50 CD34")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank) %>% 
  bind_rows(hist.plotTib)

p5 <- hist.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 CD34", "Top 50 CD34", "Next 50 CD14", "Top 50 CD14"))) %>% 
  ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
  geom_density_ridges(alpha = 0.8, quantile_lines = T, quantiles=2) +
  scale_fill_manual(values = c(c("#ADBEFF", "#3358FF", "#F9D3B4", "#EE8833"))) +
  theme(legend.position = "")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))

ggsave(p5, filename = paste0(outputFolder,"/gonzalezProportionDensity_lowComplexity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
