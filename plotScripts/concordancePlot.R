theme_set(theme_classic())

##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib   <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))
  diffpeaks <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.tsv'))
  outputFolder <- here('plots', 'concordance')
} else {
  bigtib    <- read_tsv(cmdargs[1])
  diffpeaks <- read_tsv(cmdargs[2])
  outputLFolder  <- cmdargs[3]
}

set.seed(2021)

##massage and prepare data from bigtib

###################
## RA high dose ##

#pull out the relevant values 
RA.bigtib <- bigtib %>% 
  dplyr::select(gene_name, ensg, peak, `RA-high_isDeGene`,
                `EtOH-nlDensity_avgTPM`, `RA-high_avgTPM`, `EtOH-nlDensity-avgNormFragmentCounts`, `RA-high-avgNormFragmentCounts`,
                `RA-med-avgNormFragmentCounts`, `RA-low-avgNormFragmentCounts`)

#calculate log2FC for tpm values between RA/etOH
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`RA-high_avgTPM` = `RA-high_avgTPM` + 1) %>% 
  dplyr::mutate(log2FC = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))


#annotate max peak value between RA/etOH 
RA.bigtib <- RA.bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(RA.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`RA-high-avgNormFragmentCounts`))

#annotate peak change direction
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(peakDir = case_when(`RA-high-avgNormFragmentCounts` > `EtOH-nlDensity-avgNormFragmentCounts` ~ "UP",
                                    TRUE ~ "DOWN"))

##get the differential peak identifiers for RA-high
RA.diffpeaks <- diffpeaks %>% 
  dplyr::mutate(peak = paste(chrom,startLocs,endLocs,sep=":")) %>% 
  dplyr::select(peak, `RA-high-isDiffPeak`) %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE)

##annotate if peaks in bigtib are diffpeaks
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(isDiffPeak = peak %in% RA.diffpeaks$peak)

#filter out max peak val across rows < 30
RA.bigtib <- RA.bigtib %>%
  dplyr::filter(RA.maxVal > 30)

#count number of peak combinations
RA.bigtib <- RA.bigtib %>% 
  group_by(gene_name) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & isDiffPeak == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & isDiffPeak == FALSE))

#get counts of all peaks per gene and sig peaks
plotTib <- RA.bigtib %>% 
  dplyr::select(gene_name, ensg, log2FC, Up_Sig, Up_NotSig, Down_Sig, Down_NotSig) %>% unique()

#change to proportions
plotTib <- plotTib %>% 
  rowwise() %>% 
  dplyr::mutate(totalPeaks = sum(Up_Sig, Up_NotSig, Down_Sig, Down_NotSig)) %>% 

  dplyr::mutate(Up_Sig = Up_Sig/totalPeaks) %>% 
  dplyr::mutate(Up_NotSig = Up_NotSig/totalPeaks) %>% 
  dplyr::mutate(Down_Sig = Down_Sig/totalPeaks) %>% 
  dplyr::mutate(Down_NotSig = Down_NotSig/totalPeaks) %>% 
  
  dplyr::select(- totalPeaks)


#get counts for number of peaks per gene
peakCounts <- RA.bigtib %>% dplyr::count(gene_name)

#annotate gene complexity type based on defintions from gonzalez et al, 2015
plotTib <- left_join(plotTib, peakCounts) %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::mutate(lociComplexity = case_when(numPeaks > 7 ~ "high",
                                           TRUE ~ "low"))


#make table into long format using dplyr gather
plotTib <- plotTib %>% gather("peakType", "proportion", -gene_name, -ensg, -log2FC, -numPeaks, -lociComplexity)

#change plotTib to facet 
plotTib <- plotTib %>% 
  tidyr::separate(peakType, into = c("direction", "significantPeak"))

#for clarity, filter out non significant peaks

Sig.plotTib <- plotTib %>% 
  dplyr::filter(significantPeak == "Sig")

DEGs <- RA.bigtib %>% 
  dplyr::filter(`RA-high_isDeGene` == 1)
DEGs <- DEGs$gene_name %>% unique()

hi.up.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter((lociComplexity == "high" & direction == "Up")) %>% 
  dplyr::mutate(DEG = gene_name %in% DEGs) %>% 
  dplyr::mutate(anno = case_when(log2FC > 0 & DEG == TRUE & proportion > 0.0 ~ "concordant",
                                 log2FC > 0 & DEG == TRUE & proportion == 0.0 ~ "not concordant",
                                 TRUE ~ "other"))
  
p1 <- ggscatter(hi.up.Sig.plotTib, x = "log2FC", y = "proportion", color = "anno",
                palette = c("#FE9F6D", "#8C2981", "#000000"), size = 1.2) +
  geom_hline(yintercept = 0.008) + geom_vline(xintercept = 0.5) + theme(legend.position = "")

outputLoc <- paste0(outputFolder,"/RAupScatter_0cutoff.svg")
ggsave(plot = p1, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

up.concordant.genes <- hi.up.Sig.plotTib %>% 
  dplyr::filter(anno == "concordant")
up.concordant.genes <- up.concordant.genes$gene_name

up.notconcordant.genes <- hi.up.Sig.plotTib %>% 
  dplyr::filter(anno == "not concordant")
up.notconcordant.genes <- up.notconcordant.genes$gene_name

hi.down.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter((lociComplexity == "high" & direction == "Down")) %>% 
  dplyr::mutate(DEG = gene_name %in% DEGs) %>% 
  dplyr::mutate(anno = case_when(log2FC < 0 & DEG == TRUE & proportion > 0.0 ~ "concordant",
                                 log2FC < 0 & DEG == TRUE & proportion == 0.0 ~ "not concordant",
                                 TRUE ~ "other"))

p2 <- ggscatter(hi.down.Sig.plotTib, x = "log2FC", y = "proportion", color = "anno",
                palette = c("#DE4968", "#3B0f70", "#000000"), size = 1.2) +
  geom_hline(yintercept = 0.0058) + geom_vline(xintercept = -0.5) + theme(legend.position = "") +
  ylim(0,1)

outputLoc <- paste0(outputFolder,"/RAdownScatter_0cutoff.svg")
ggsave(plot = p2, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

down.concordant.genes <- hi.down.Sig.plotTib %>% 
  dplyr::filter(anno == "concordant")
down.concordant.genes <- down.concordant.genes$gene_name

down.notconcordant.genes <- hi.down.Sig.plotTib %>% 
  dplyr::filter(anno == "not concordant")
down.notconcordant.genes <- down.notconcordant.genes$gene_name

##
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(anno = case_when(gene_name %in% up.concordant.genes ~ "Up concordant",
                                 gene_name %in% up.notconcordant.genes ~ "Up not concordant",
                                 gene_name %in% down.concordant.genes ~ "Down concordant",
                                 gene_name %in% down.notconcordant.genes ~ "Down not concordant",
                                 TRUE ~ "other")) %>% 
  dplyr::mutate(anno = factor(anno, levels = c("Down concordant",
                                               "Down not concordant",
                                               "Up concordant",
                                               "Up not concordant",
                                               "other")))

p3 <-  RA.bigtib %>% 
  ggplot(aes(x=anno, y=`EtOH-nlDensity-avgNormFragmentCounts`, color = anno)) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width=0, justification = -0.3, point_colour=NA) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(size = 0.2, alpha = 0.4, position = position_jitter(seed=1, width = 0.1)) + 
  ylim(0,200) +
  scale_color_manual(values = c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000")) + 
  scale_fill_manual(values =c("#DE4968", "#3B0f70","#FE9F6D", "#8C2981", "#000000")) + 
  theme(legend.position = "")

outputLoc <- paste0(outputFolder,"/RA_etOHrainCloud_0cutoff.svg")
ggsave(plot = p3, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

###perform stats (to be added in illustrator)###
down.not.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Down not concordant")
down.not.concordant <- down.not.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

down.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Down concordant")
down.concordant <- down.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

ks.test(down.not.concordant, down.concordant)

up.not.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Up not concordant")
up.not.concordant <- up.not.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

up.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Up concordant")
up.concordant <- up.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

ks.test(up.not.concordant, up.concordant)

p4 <- RA.bigtib %>% 
  ggplot(aes(x=anno, y=`RA-high-avgNormFragmentCounts`, color = anno)) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width=0, justification = -0.3, point_colour=NA) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(size = 0.5, alpha = 0.3, position = position_jitter(seed=1, width = 0.1)) + 
  ylim(0,200) +
  scale_color_manual(values = c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000")) + 
  scale_fill_manual(values =c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000")) + 
  theme(legend.position = "")

outputLoc <- paste0(outputFolder,"/RA_rarainCloud_0cutoff.svg")
ggsave(plot = p4, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

down.not.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Down not concordant")
down.not.concordant <- down.not.concordant$`RA-high-avgNormFragmentCounts`

down.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Down concordant")
down.concordant <- down.concordant$`RA-high-avgNormFragmentCounts`

ks.test(down.not.concordant, down.concordant)

up.not.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Up not concordant")
up.not.concordant <- up.not.concordant$`RA-high-avgNormFragmentCounts`

up.concordant <- RA.bigtib %>% 
  dplyr::filter(anno == "Up concordant")
up.concordant <- up.concordant$`RA-high-avgNormFragmentCounts`

ks.test(up.not.concordant, up.concordant)

#check the relative peak complexity between groups
peakComplexity.up.tib <- hi.up.Sig.plotTib %>% 
  dplyr::select(anno, numPeaks) %>% 
  dplyr::mutate(anno = case_when(anno == "concordant" ~ "Up concordant",
                                 anno == "not concordant" ~ "Up not concordant",
                                 TRUE ~ "other"))

peakComplexity.down.tib <-  hi.down.Sig.plotTib %>% 
  dplyr::select(anno, numPeaks) %>% 
  dplyr::mutate(anno = case_when(anno == "concordant" ~ "Down concordant",
                                 anno == "not concordant" ~ "Down not concordant",
                                 TRUE ~ "other"))

peakComplexity.tib <- bind_rows(peakComplexity.up.tib, peakComplexity.down.tib) %>% 
  dplyr::mutate(anno = factor(anno, levels = c("Down concordant",
                                               "Down not concordant",
                                               "Up concordant",
                                               "Up not concordant",
                                               "other")))

p5 <- ggboxplot(peakComplexity.tib, x= "anno", y = "numPeaks", color = "anno",
                palette = c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000"), outlier.shape = NA) + ylim(0,40)

outputLoc <- paste0(outputFolder, "/RAcomplexitydist_0cutoff.svg")
ggsave(plot = p5, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

###Stats!! ###
up.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Up concordant")

up.concordant <- up.concordant$numPeaks

up.not.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Up not concordant")

up.not.concordant <- up.not.concordant$numPeaks

down.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Down concordant")

down.concordant <- down.concordant$numPeaks

down.not.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Down not concordant")

down.not.concordant <- down.not.concordant$numPeaks


ks.test(up.not.concordant, down.not.concordant)
ks.test(up.concordant, down.concordant)
ks.test(up.not.concordant, up.concordant)
ks.test(down.concordant, down.not.concordant)

#perform some stats on distributions
other <- peakComplexity.tib %>% 
  dplyr::filter(anno == "other")
other <- other$numPeaks

up.not.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Up not concordant")
up.not.concordant <- up.not.concordant$numPeaks

up.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Up concordant")
up.concordant <- up.concordant$numPeaks

down.not.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Down not concordant")
down.not.concordant <- down.not.concordant$numPeaks

down.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Down concordant")
down.concordant <- down.concordant$numPeaks

ks.test(up.concordant, down.concordant)
#####################
## TGFb high dose ##

#pull out the relevant values 
TGFb.bigtib <- bigtib %>% 
  dplyr::select(gene_name, peak, `TGFb-high_isDeGene`,
                `EtOH-nlDensity_avgTPM`, `TGFb-high_avgTPM`, `EtOH-nlDensity-avgNormFragmentCounts`, `TGFb-high-avgNormFragmentCounts`,
                `TGFb-med-avgNormFragmentCounts`, `TGFb-low-avgNormFragmentCounts`)

#calculate log2FC for tpm values between TGFb/etOH
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`TGFb-high_avgTPM` = `TGFb-high_avgTPM` + 1) %>% 
  dplyr::mutate(log2FC = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))


#annotate max peak value between TGFb/etOH 
TGFb.bigtib <- TGFb.bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(TGFb.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`TGFb-high-avgNormFragmentCounts`))

#annotate peak change direction
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(peakDir = case_when(`TGFb-high-avgNormFragmentCounts` > `EtOH-nlDensity-avgNormFragmentCounts` ~ "UP",
                                    TRUE ~ "DOWN"))

##get the differential peak identifiers for TGFb-high
TGFb.diffpeaks <- diffpeaks %>% 
  dplyr::mutate(peak = paste(chrom,startLocs,endLocs,sep=":")) %>% 
  dplyr::select(peak, `TGFb-high-isDiffPeak`) %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE)

##annotate if peaks in bigtib are diffpeaks
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffpeaks$peak)

#filter out max peak val across rows < 30
TGFb.bigtib <- TGFb.bigtib %>%
  dplyr::filter(TGFb.maxVal > 30)

#count number of peak combinations
TGFb.bigtib <- TGFb.bigtib %>% 
  group_by(gene_name) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & isDiffPeak == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & isDiffPeak == FALSE))

#get counts of all peaks per gene and sig peaks
plotTib <- TGFb.bigtib %>% 
  dplyr::select(gene_name, log2FC, Up_Sig, Up_NotSig, Down_Sig, Down_NotSig) %>% unique()

#change to proportions
plotTib <- plotTib %>% 
  dplyr::mutate(totalPeaks = sum(Up_Sig, Up_NotSig, Down_Sig, Down_NotSig)) %>% 
  
  dplyr::mutate(Up_Sig = Up_Sig/totalPeaks) %>% 
  dplyr::mutate(Up_NotSig = Up_NotSig/totalPeaks) %>% 
  dplyr::mutate(Down_Sig = Down_Sig/totalPeaks) %>% 
  dplyr::mutate(Down_NotSig = Down_NotSig/totalPeaks) %>% 
  
  dplyr::select(- totalPeaks)


#get counts for number of peaks per gene
peakCounts <- TGFb.bigtib %>% dplyr::count(gene_name)

#annotate gene complexity type based on defintions from gonzalez et al, 2015
plotTib <- left_join(plotTib, peakCounts) %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::mutate(lociComplexity = case_when(numPeaks > 7 ~ "high",
                                           TRUE ~ "low"))


#make table into long format using dplyr gather
plotTib <- plotTib %>% gather("peakType", "proportion", -gene_name, -log2FC, -numPeaks, -lociComplexity)

#change plotTib to facet 
plotTib <- plotTib %>% 
  tidyr::separate(peakType, into = c("direction", "significantPeak"))

#for clarity, filter out non significant peaks

Sig.plotTib <- plotTib %>% 
  dplyr::filter(significantPeak == "Sig")

DEGs <- TGFb.bigtib %>% 
  dplyr::filter(`TGFb-high_isDeGene` == 1)
DEGs <- DEGs$gene_name %>% unique()

hi.up.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter((lociComplexity == "high" & direction == "Up")) %>% 
  dplyr::mutate(DEG = gene_name %in% DEGs) %>% 
  dplyr::mutate(anno = case_when(log2FC > 0 & DEG == TRUE & proportion > 0.0 ~ "concordant",
                                 log2FC > 0 & DEG == TRUE & proportion == 0.0 ~ "not concordant",
                                 TRUE ~ "other"))

p6 <- ggscatter(hi.up.Sig.plotTib, x = "log2FC", y = "proportion", color = "anno",
                palette = c("#FE9F6D", "#8C2981", "#000000"), size = 1.2) +
  geom_hline(yintercept = 0.008) + geom_vline(xintercept = 0.5) + theme(legend.position = "") +
  ylim(0,1)

outputLoc <- paste0(outputFolder,"/TGFbupScatter_0cutoff.svg")
ggsave(plot = p6, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

up.concordant.genes <- hi.up.Sig.plotTib %>% 
  dplyr::filter(anno == "concordant")
up.concordant.genes <- up.concordant.genes$gene_name

up.notconcordant.genes <- hi.up.Sig.plotTib %>% 
  dplyr::filter(anno == "not concordant")
up.notconcordant.genes <- up.notconcordant.genes$gene_name

hi.down.Sig.plotTib <- Sig.plotTib %>% 
  dplyr::filter((lociComplexity == "high" & direction == "Down")) %>% 
  dplyr::mutate(DEG = gene_name %in% DEGs) %>% 
  dplyr::mutate(anno = case_when(log2FC < 0 & DEG == TRUE & proportion > 0.0 ~ "concordant",
                                 log2FC < 0 & DEG == TRUE & proportion == 0.0 ~ "not concordant",
                                 TRUE ~ "other"))

p7 <- ggscatter(hi.down.Sig.plotTib, x = "log2FC", y = "proportion", color = "anno",
                palette = c("#DE4968", "#3B0f70", "#000000"), size = 1.2) +
  geom_hline(yintercept = 0.008) + geom_vline(xintercept = -0.5) + theme(legend.position = "") +
  ylim(0,1)

outputLoc <- paste0(outputFolder,"/TGFbdownScatter_0cutoff.svg")
ggsave(plot = p7, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

down.concordant.genes <- hi.down.Sig.plotTib %>% 
  dplyr::filter(anno == "concordant")
down.concordant.genes <- down.concordant.genes$gene_name

down.notconcordant.genes <- hi.down.Sig.plotTib %>% 
  dplyr::filter(anno == "not concordant")
down.notconcordant.genes <- down.notconcordant.genes$gene_name

##
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(anno = case_when(gene_name %in% up.concordant.genes ~ "Up concordant",
                                 gene_name %in% up.notconcordant.genes ~ "Up not concordant",
                                 gene_name %in% down.concordant.genes ~ "Down concordant",
                                 gene_name %in% down.notconcordant.genes ~ "Down not concordant",
                                 TRUE ~ "other")) %>% 
  dplyr::mutate(anno = factor(anno, levels = c("Down concordant",
                                               "Down not concordant",
                                               "Up concordant",
                                               "Up not concordant",
                                               "other")))

p8 <-  TGFb.bigtib %>% 
  ggplot(aes(x=anno, y=`EtOH-nlDensity-avgNormFragmentCounts`, color = anno, fill = anno)) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width=0, justification = -0.3, point_colour=NA) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(size = 0.2, alpha = 0.4, position = position_jitter(seed=1, width = 0.1)) + 
  ylim(0,200) +
  scale_color_manual(values = c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000")) + 
  scale_fill_manual(values =c("#DE4968", "#3B0f70","#FE9F6D", "#8C2981", "#000000")) + 
  theme(legend.position = "")

outputLoc <- paste0(outputFolder,"/TGFb_etOHrainCloud_0cutoff.svg")
ggsave(plot = p8, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

###perform stats (to be added in illustrator)###
down.not.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Down not concordant")
down.not.concordant <- down.not.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

down.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Down concordant")
down.concordant <- down.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

ks.test(down.not.concordant, down.concordant)

up.not.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Up not concordant")
up.not.concordant <- up.not.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

up.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Up concordant")
up.concordant <- up.concordant$`EtOH-nlDensity-avgNormFragmentCounts`

ks.test(up.not.concordant, up.concordant)

p9 <- TGFb.bigtib %>% 
  ggplot(aes(x=anno, y=`TGFb-high-avgNormFragmentCounts`, color = anno, fill = anno)) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width=0, justification = -0.3, point_colour=NA) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(size = 0.5, alpha = 0.3, position = position_jitter(seed=1, width = 0.1)) + 
  ylim(0,200) +
  scale_color_manual(values = c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000")) + 
  scale_fill_manual(values =c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000")) + 
  theme(legend.position = "")

outputLoc <- paste0(outputFolder,"/TGFb_tgfbrainCloud_0cutoff.svg")
ggsave(plot = p9, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

###stats
down.not.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Down not concordant")
down.not.concordant <- down.not.concordant$`TGFb-high-avgNormFragmentCounts`

down.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Down concordant")
down.concordant <- down.concordant$`TGFb-high-avgNormFragmentCounts`

ks.test(down.not.concordant, down.concordant)

up.not.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Up not concordant")
up.not.concordant <- up.not.concordant$`TGFb-high-avgNormFragmentCounts`

up.concordant <- TGFb.bigtib %>% 
  dplyr::filter(anno == "Up concordant")
up.concordant <- up.concordant$`TGFb-high-avgNormFragmentCounts`

ks.test(up.not.concordant, up.concordant)

#check the relative peak complexity between groups
peakComplexity.up.tib <- hi.up.Sig.plotTib %>% 
  dplyr::select(anno, numPeaks) %>% 
  dplyr::mutate(anno = case_when(anno == "concordant" ~ "Up concordant",
                                 anno == "not concordant" ~ "Up not concordant",
                                 TRUE ~ "other"))

peakComplexity.down.tib <-  hi.down.Sig.plotTib %>% 
  dplyr::select(anno, numPeaks) %>% 
  dplyr::mutate(anno = case_when(anno == "concordant" ~ "Down concordant",
                                 anno == "not concordant" ~ "Down not concordant",
                                 TRUE ~ "other"))

peakComplexity.tib <- bind_rows(peakComplexity.up.tib, peakComplexity.down.tib) %>% 
  dplyr::mutate(anno = factor(anno, levels = c("Down concordant",
                                               "Down not concordant",
                                               "Up concordant",
                                               "Up not concordant",
                                               "other")))

p10 <- ggboxplot(peakComplexity.tib, x= "anno", y = "numPeaks", color = "anno",
                palette = c("#DE4968", "#3B0f70", "#FE9F6D", "#8C2981", "#000000"), outlier.shape = NA) + ylim(0,40)

outputLoc <- paste0(outputFolder, "/TGFbcomplexitydist_0cutoff.svg")
ggsave(plot = p10, filename = outputLoc, width=5, height = 5, units = "in", dpi = 300)

###Stats!! ###
up.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Up concordant")

up.concordant <- up.concordant$numPeaks

up.not.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Up not concordant")

up.not.concordant <- up.not.concordant$numPeaks

down.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Down concordant")

down.concordant <- down.concordant$numPeaks

down.not.concordant <- peakComplexity.tib %>% 
  dplyr::filter(anno == "Down not concordant")

down.not.concordant <- down.not.concordant$numPeaks


ks.test(up.not.concordant, down.not.concordant)
ks.test(up.concordant, down.concordant)
ks.test(up.not.concordant, up.concordant)
ks.test(down.concordant, down.not.concordant)
