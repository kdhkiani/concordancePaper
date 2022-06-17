theme_set(theme_classic())

sanfordMakeProportionPlot <- function(bigtib_path, diffpeaks_path, outputLoc){
  
  bigtib <- read_tsv(bigtib_path)
  diffpeaks <- read_tsv(diffpeaks_path)
  
  ##massage and prepare data from bigtib
  
  ###################
  ## RA high dose ##
  
  #pull out the relevant values 
  RA.bigtib <- bigtib %>% 
    dplyr::select(gene_name, peak, `RA-high_isDeGene`,
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
  
  #filter out max peak val across lows < 20
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
  peakCounts <- RA.bigtib %>% dplyr::count(gene_name)
  
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
  
  ### get rank values
  
  #now get ranks
  top50 <- Sig.plotTib %>% 
    arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  top50 <- top50$gene_name
  
  top100 <- Sig.plotTib %>% 
    arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  top100 <- top100$gene_name
  
  
  bot50 <- Sig.plotTib %>% 
    arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  bot50 <- bot50$gene_name
  
  bot100<- Sig.plotTib %>% 
    arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  bot100 <- bot100$gene_name
  
  ##annotate based on rank position
  Sig.plotTib <- Sig.plotTib %>% 
    dplyr::mutate(geneRank = case_when(gene_name %in% top50 ~ "50 most upregulated",
                                       gene_name %in% top100 ~ "Next 50 most upregulated",
                                       gene_name %in% bot50 ~ "50 most downregulated",
                                       gene_name %in% bot100 ~ "Next 50 most downregulated",
                                       TRUE ~ "other"))
  
  hi.Sig.plotTib <- Sig.plotTib %>% 
    dplyr::filter(lociComplexity == "high")
  
  ###plot
  p1 <- ggscatter(hi.Sig.plotTib, x = "log2FC", y = "proportion",
                  facet.by = c("direction"), color = "geneRank",
                  palette = c("#808782", "#C1292E", "#CACECB", "#E6898C", "#000000"),
                  alpha = 0.8, size = 1.5)
  
  ggsave(p1, filename = paste0(outputLoc,"RAProportionPlot_highComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)
  
  lo.Sig.plotTib <- Sig.plotTib %>% 
    dplyr::filter(lociComplexity == "low")
  
  ###plot
  p2 <- ggscatter(lo.Sig.plotTib, x = "log2FC", y = "proportion",
                  facet.by = c("direction"), color = "geneRank",
                  palette = c("#808782", "#C1292E", "#CACECB", "#E6898C", "#000000"),
                  alpha = 0.8, size = 1.5)
  
  ggsave(p2, filename = paste0(outputLoc,"RAProportionPlot_lowComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)
  
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Up" & lociComplexity == "high") %>% 
    dplyr::filter(geneRank == "50 most upregulated" | geneRank == "Next 50 most upregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank)
  
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Down" & lociComplexity == "high") %>% 
    dplyr::filter(geneRank == "50 most downregulated" | geneRank == "Next 50 most downregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank) %>% 
    bind_rows(hist.plotTib)
  
  p3 <- hist.plotTib %>% 
    dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 most downregulated", "50 most downregulated", "Next 50 most upregulated", "50 most upregulated"))) %>% 
    ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
    geom_density_ridges(alpha = 0.8, quantile_lines=TRUE, quantiles=2) +
    scale_fill_manual(values = c(c("#CACECB", "#808782","#E6898C", "#C1292E"))) +
    theme(legend.position = "") +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))
  
  ggsave(p3, filename = paste0(outputLoc,"raProportionDensity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
  
  ##
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Up" & lociComplexity == "low") %>% 
    dplyr::filter(geneRank == "50 most upregulated" | geneRank == "Next 50 most upregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank)
  
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Down" & lociComplexity == "low") %>% 
    dplyr::filter(geneRank == "50 most downregulated" | geneRank == "Next 50 most downregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank) %>% 
    bind_rows(hist.plotTib)
  
  p4 <- hist.plotTib %>% 
    dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 most downregulated", "50 most downregulated", "Next 50 most upregulated", "50 most upregulated"))) %>% 
    ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
    geom_density_ridges(alpha = 0.8, quantile_lines=TRUE, quantiles=2) +
    scale_fill_manual(values = c(c("#CACECB", "#808782","#E6898C", "#C1292E"))) +
    theme(legend.position = "") +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))
  
  ggsave(p4, filename = paste0(outputLoc,"RAProportionDensity_lowComplexity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
  
  
  ##########################
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
  
  #filter out max peak val across lows < 20
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
  
  ### get rank values
  
  #now get ranks
  top50 <- Sig.plotTib %>% 
    arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  top50 <- top50$gene_name
  
  top100 <- Sig.plotTib %>% 
    arrange(desc(log2FC)) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  top100 <- top100$gene_name
  
  
  bot50 <- Sig.plotTib %>% 
    arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(1:100) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  bot50 <- bot50$gene_name
  
  bot100<- Sig.plotTib %>% 
    arrange(log2FC) %>% group_by(lociComplexity) %>% dplyr::slice(101:200) %>% 
    dplyr::select(gene_name, log2FC, lociComplexity) %>% unique()
  bot100 <- bot100$gene_name
  
  ##annotate based on rank position
  Sig.plotTib <- Sig.plotTib %>% 
    dplyr::mutate(geneRank = case_when(gene_name %in% top50 ~ "50 most upregulated",
                                       gene_name %in% top100 ~ "Next 50 most upregulated",
                                       gene_name %in% bot50 ~ "50 most downregulated",
                                       gene_name %in% bot100 ~ "Next 50 most downregulated",
                                       TRUE ~ "other"))
  
  hi.Sig.plotTib <- Sig.plotTib %>% 
    dplyr::filter(lociComplexity == "high")
  
  ###plot
  p5 <- ggscatter(hi.Sig.plotTib, x = "log2FC", y = "proportion",
                  facet.by = c("direction"), color = "geneRank",
                  palette = c("#808782", "#091e5f", "#CACECB", "#113cbb", "#000000"),
                  alpha = 0.8, size = 1.5) + ylim(0,1)
  
  ggsave(p5, filename = paste0(outputLoc,"TGFbProportionPlot_highComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)
  
  lo.Sig.plotTib <- Sig.plotTib %>% 
    dplyr::filter(lociComplexity == "low")
  
  ###plot
  p6 <- ggscatter(lo.Sig.plotTib, x = "log2FC", y = "proportion",
                  facet.by = c("direction"), color = "geneRank",
                  palette = c("#808782", "#091e5f", "#CACECB", "#113cbb", "#000000"),
                  alpha = 0.8, size = 1.5)
  
  ggsave(p6, filename = paste0(outputLoc,"TGFbProportionPlot_lowComplexity.svg"),  width = 7.5, height = 3.75, units = "in", dpi = 300)
  
  
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Up" & lociComplexity == "high") %>% 
    dplyr::filter(geneRank == "50 most upregulated" | geneRank == "Next 50 most upregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank)
  
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Down" & lociComplexity == "high") %>% 
    dplyr::filter(geneRank == "50 most downregulated" | geneRank == "Next 50 most downregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank) %>% 
    bind_rows(hist.plotTib)
  
  p7 <- hist.plotTib %>% 
    dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 most downregulated", "50 most downregulated", "Next 50 most upregulated", "50 most upregulated"))) %>% 
    ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
    geom_density_ridges(alpha = 0.8, quantile_lines=TRUE, quantiles=2) +
    scale_fill_manual(values = c("#CACECB", "#808782", "#113cbb", "#091e5f")) +
    theme(legend.position = "") +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))
  
  ggsave(p7, filename = paste0(outputLoc,"tgfbProportionDensity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
  
  ##
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Up" & lociComplexity == "low") %>% 
    dplyr::filter(geneRank == "50 most upregulated" | geneRank == "Next 50 most upregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank)
  
  tmp <- Sig.plotTib %>% 
    dplyr::filter(direction == "Down" & lociComplexity == "low") %>% 
    dplyr::filter(geneRank == "50 most downregulated" | geneRank == "Next 50 most downregulated")
  
  hist.plotTib <- tmp %>% 
    dplyr::select(proportion, geneRank) %>% 
    bind_rows(hist.plotTib)
  
  p8 <- hist.plotTib %>% 
    dplyr::mutate(geneRank = factor(geneRank, levels = c("Next 50 most downregulated", "50 most downregulated", "Next 50 most upregulated", "50 most upregulated"))) %>% 
    ggplot( aes(x= proportion, y = geneRank, fill = geneRank)) +
    geom_density_ridges(alpha = 0.8, quantile_lines=TRUE, quantiles=2) +
    scale_fill_manual(values = c("#CACECB", "#808782", "#113cbb", "#091e5f")) +
    theme(legend.position = "") +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.1), limits = c(0,1))
  
  ggsave(p8, filename = paste0(outputLoc,"TGFbProportionDensity_lowComplexity_median.svg"),  width = 5, height = 5, units = "in", dpi = 300)
  

  
}
