theme_set(theme_classic())

########
#using source data from paper website
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib       <- read_tsv(here('extractedData', 'Ramirez', 'ramirez_genePeak_1to1mapping_mergedDist50.tsv'))
  outputFolder <- here('plots', 'proportion')
} else {
  bigtib         <- cmdargs[1]
  outputFolder   <- cmdargs[2]
}

##massage and prepare data from bigtib

#drop any peaks with no genes assigned
bigtib <- drop_na(bigtib)

#find average across replicates
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(HL60_avgTPM = mean(c(SRR3214253_HL60_tpm, SRR3214254_HL60_tpm, SRR3214255_HL60_tpm))) %>% 
  dplyr::mutate(monocyte_avgTPM = mean(c(SRR3214325_monocyte_tpm, SRR3214326_monocyte_tpm, SRR3214327_monocyte_tpm))) %>% 
  dplyr::mutate(neutrophil_12hr_avgTPM = mean(c(SRR3214289_neutrophil12hr_tpm, SRR3214290_neutrophil12hr_tpm, SRR3214291_neutrophil12hr_tpm))) %>% 
  dplyr::mutate(neutrophil_120hr_avgTPM = mean(c(SRR3214301_neutrophil120hr_tpm, SRR3214302_neutrophil120hr_tpm, SRR3214303_neutrophil120hr_tpm))) %>% 
  dplyr::mutate(monomac_avgTPM = mean(c(SRR3214343_monomac_tpm, SRR3214344_monomac_tpm, SRR3214345_monomac_tpm))) %>% 
  
  
  dplyr::mutate(HL60_avgNormFragmentCounts = mean(c(SRR3213963_HL60_normFragmentCounts, SRR3213964_HL60_normFragmentCounts, SRR3213965_HL60_normFragmentCounts))) %>% 
  dplyr::mutate(monocyte_avgNormFragmentCounts = mean(c(SRR3214035_monocyte_normFragmentCounts, SRR3214036_monocyte_normFragmentCounts, SRR3214037_monocyte_normFragmentCounts))) %>% 
  dplyr::mutate(neutrophil_12hr_avgNormFragmentCounts = mean(c(SRR3213999_neutrophil_12hr_normFragmentCounts, SRR3214000_neutrophil_12hr_normFragmentCounts, SRR3214001_neutrophil_12hr_normFragmentCounts))) %>% 
  dplyr::mutate(neutrophil_120hr_avgNormFragmentCounts = mean(c(SRR3214011_neutrophil_120hr_normFragmentCounts, SRR3214012_neutrophil_120hr_normFragmentCounts, SRR3214013_neutrophil_120hr_normFragmentCounts))) %>% 
  dplyr::mutate(monomac_avgNormFragmentCounts = mean(c(SRR3214053_monomac_168hr_normFragmentCounts, SRR3214054_monomac_168hr_normFragmentCounts, SRR3214055_monomac_168hr_normFragmentCounts)))

#annotate max peak value between HL60 and monocyte conditions
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(maxValATAC = max(HL60_avgNormFragmentCounts, monocyte_avgNormFragmentCounts, neutrophil_12hr_avgNormFragmentCounts, neutrophil_120hr_avgNormFragmentCounts, monomac_avgNormFragmentCounts)) %>% 
  dplyr::mutate(maxValRNA = max(HL60_avgTPM, monocyte_avgTPM, neutrophil_12hr_avgTPM, neutrophil_120hr_avgTPM, monomac_avgTPM))

#create histogram of distribution of max value between two samples of RNA
p1 <- ggplot(bigtib, aes(x = maxValRNA)) + 
  geom_histogram(fill = "light gray", color = "black", alpha = 0.8, binwidth = 1) + 
  xlim(-1,50) + geom_vline(xintercept = 5) 

#add 1 to avoid infinite values when doing log2 fold changes
bigtib <- bigtib %>% 
  dplyr::mutate(HL60_avgTPM = HL60_avgTPM + 1) %>% 
  dplyr::mutate(monocyte_avgTPM = monocyte_avgTPM + 1) %>% 
  dplyr::mutate(neutrophil_12hr_avgTPM = neutrophil_12hr_avgTPM + 1) %>% 
  dplyr::mutate(neutrophil_120hr_avgTPM = neutrophil_120hr_avgTPM + 1) %>% 
  dplyr::mutate(monomac_avgTPM = monomac_avgTPM + 1) %>% 
  
  dplyr::mutate(HL60_avgNormFragmentCounts = HL60_avgNormFragmentCounts + 1) %>% 
  dplyr::mutate(monocyte_avgNormFragmentCounts = monocyte_avgNormFragmentCounts + 1) %>% 
  dplyr::mutate(neutrophil_12hr_avgNormFragmentCounts = neutrophil_12hr_avgNormFragmentCounts + 1) %>% 
  dplyr::mutate(neutrophil_120hr_avgNormFragmentCounts = neutrophil_120hr_avgNormFragmentCounts + 1) %>% 
  dplyr::mutate(monomac_avgNormFragmentCounts = monomac_avgNormFragmentCounts + 1)

#create histogram of distribution of max value between two samples of RNA
p2 <- ggplot(bigtib, aes(x = maxValATAC)) + 
  geom_histogram(fill = "light gray", color = "black", alpha = 0.8, binwidth = 1) + 
  xlim(-1,50) + geom_vline(xintercept = 5) 

#outputLoc <- paste0(outputFolder, "/maxRNAhist.svg")
#ggsave(plot = p1, filename = outputLoc, height = 5, width = 5, units = "in", dpi = 300)

#remove genes where max rna between cell populations is less than 5 tpm
bigtib <- bigtib %>% 
  dplyr::filter(maxValRNA > 5) %>% 
  dplyr::filter(maxValATAC > 5)

###First, look at difference between monocyte-derived macrophages and HL60 cells###


#calculate log2FC for accessibility and create a filtering column
bigtib <- bigtib %>% 
  dplyr::mutate(log2FCatac = log2(monomac_avgNormFragmentCounts/HL60_avgNormFragmentCounts)) %>% 
  dplyr::mutate(peakCutoff = abs(log2FCatac) > 0.58)

#annotate peak change direction
bigtib <- bigtib %>% 
  dplyr::mutate(peakDir = case_when(monomac_avgNormFragmentCounts > HL60_avgNormFragmentCounts ~ "UP",
                                    TRUE ~ "DOWN"))
#determine log2FC for expression
bigtib <- bigtib %>% 
  dplyr::mutate(log2FC = log2(monomac_avgTPM/HL60_avgTPM))

#arrange based on log2FC expression values
bigtib <- bigtib %>% 
  dplyr::arrange(desc(log2FC))

#count number of peak combinations
bigtib <- bigtib %>% 
  group_by(ensg) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & peakCutoff == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & peakCutoff == FALSE))

#get counts for number of peaks per gene
peakCounts <- bigtib %>% dplyr::count(ensg)

#get counts of all peaks per gene and sig peaks
plotTib <- bigtib %>% 
  dplyr::select(ensg, log2FC, Up_Sig, Up_NotSig, Down_Sig, Down_NotSig) %>% unique()

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

plotTib <- plotTib %>% gather("peakType", "proportion", -ensg, -log2FC, -numPeaks, -lociComplexity)


#change plotTib to attempt facetting 
plotTib <- plotTib %>% 
  tidyr::separate(peakType, into = c("direction", "significantPeak"))


#for clarity, filter out non significant peaks

Sig.plotTib <- plotTib %>% 
  dplyr::filter(significantPeak == "Sig") 

### get rank values
top50 <- Sig.plotTib %>%
  arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>%
  dplyr::select(ensg, log2FC, lociComplexity) %>% unique()
top50 <- top50$ensg

top100 <- Sig.plotTib %>%
  arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>%
  dplyr::select(ensg, log2FC, lociComplexity) %>% unique()
top100 <- top100$ensg

bot50 <- Sig.plotTib %>%
  arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>%
  dplyr::select(ensg, log2FC, lociComplexity) %>% unique()
bot50 <- bot50$ensg

bot100<- Sig.plotTib %>%
  arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>%
  dplyr::select(ensg, log2FC, lociComplexity) %>% unique()
bot100 <- bot100$ensg

##annotate based on rank position
Sig.plotTib <- Sig.plotTib %>% 
  dplyr::mutate(geneRank = case_when(ensg %in% top50 ~ "Top 50 macrophage",
                                     ensg %in% top100 ~ "Next 50 macrophage",
                                     ensg %in% bot50 ~ "Top 50 HL60",
                                     ensg %in% bot100 ~ "Next 50 HL60",
                                     TRUE ~ "other"))

Sig.plotTib <- Sig.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Top 50 HL60", "Next 50 HL60",
                                                       "other",
                                                       "Next 50 macrophage", "Top 50 macrophage"
                                                       )))
hi.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter(lociComplexity == "high")

p2 <- ggscatter(hi.Sig.plotTib, x = "log2FC", y = "proportion",
                facet.by = c("direction"), color = "geneRank",
                palette = c("#9DC4B5", "#D9E8E2", "#000000", "#D6300A", "#A72608"),
                alpha = 0.8, size = 1.5) +
  ggtitle("Proportion of ATAC peaks differenitally accessible for top 100 untreated and IL1alpha treated genes from Weiterer et al., 2020") +
  xlab("log2 Fold Change in Expression") +
  ylab("Proportion of DHS peaks differentially accessible") +
  ylim(0,1)

ggsave(p2, filename = paste0(outputFolder,"/ramirezProportionPlot_highComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)

lo.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter(lociComplexity == "low")

p3 <- ggscatter(lo.Sig.plotTib, x = "log2FC", y = "proportion",
                facet.by = c("direction"), color = "geneRank",
                palette = c("#9DC4B5", "#D9E8E2", "#000000", "#D6300A", "#A72608"),
                alpha = 0.8, size = 1.5) +
  ggtitle("Proportion of DHS peaks differenitally accessible for top 100 CD34 and CD14 genes from Gonzalez et al., 2015") +
  xlab("log2 Fold Change in Expression") +
  ylab("Proportion of DHS peaks differentially accessible") +
  ylim(0,1)

ggsave(p3, filename = paste0(outputFolder,"/ramirezProportionPlot_lowComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Up" & lociComplexity == "high") %>% 
  dplyr::filter(geneRank == "Top 50 macrophage" | geneRank == "Next 50 macrophage")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Down" & lociComplexity == "high") %>% 
  dplyr::filter(geneRank == "Top 50 HL60" | geneRank == "Next 50 HL60")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank) %>% 
  bind_rows(hist.plotTib)

p4 <- hist.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Top 50 HL60", "Next 50 HL60",
                                                       "Next 50 macrophage", "Top 50 macrophage"))) %>% 
  ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
  geom_density_ridges(alpha = 0.8, quantile_lines = T, quantiles=2) +
  scale_fill_manual(values = c(c("#9DC4B5", "#D9E8E2","#D6300A", "#A72608"))) +
  theme(legend.position = "")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))

ggsave(p4, filename = paste0(outputFolder,"/ramirezProportionDensity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)


tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Up" & lociComplexity == "low") %>% 
  dplyr::filter(geneRank == "Top 50 macrophage" | geneRank == "Next 50 macrophage")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank)

tmp <- Sig.plotTib %>% 
  dplyr::filter(direction == "Down" & lociComplexity == "low") %>% 
  dplyr::filter(geneRank == "Top 50 HL60" | geneRank == "Next 50 HL60")

hist.plotTib <- tmp %>% 
  dplyr::select(proportion, geneRank) %>% 
  bind_rows(hist.plotTib)

p5 <- hist.plotTib %>% 
  dplyr::mutate(geneRank = factor(geneRank, levels = c("Top 50 HL60", "Next 50 HL60",
                                                       "Next 50 macrophage", "Top 50 macrophage"))) %>% 
  ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
  geom_density_ridges(alpha = 0.8, quantile_lines = T, quantiles=2) +
  scale_fill_manual(values = c(c("#9DC4B5", "#D9E8E2","#D6300A", "#A72608"))) +
  theme(legend.position = "")+
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))

ggsave(p5, filename = paste0(outputFolder,"/ramirezProportionDensity_lowComplexity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
