theme_set(theme_classic())

sanfordScatterPlots <- function(bigtib_path, diffpeaks_path, outputFolder){
  set.seed(2021)
  
  bigTib <- read_tsv(bigtib_path)
  diffPeaks <- read_tsv(diffpeaks_path)
  
  if (! "peak" %in% colnames(bigTib)){
    bigTib <- bigTib %>% 
      dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":"))
  }
  
  diffPeaks <- diffPeaks %>% 
    dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":"))
  
  subsetTib <- dplyr::select(bigTib, matches("(gene_name)|(peak)|(avgNormFragmentCounts)|(_avgTPM)")) %>% 
    dplyr::select(-contains("highDensity")) %>%
    dplyr::select(-contains("halfDensity")) %>% 
    dplyr::select(-contains("-and-"))
  
  RA.diffPeaks <- diffPeaks %>% 
    dplyr::filter(`RA-high-isDiffPeak` == TRUE) %>% 
    dplyr::select(peak)
  
  TGFb.diffPeaks <- diffPeaks %>% 
    dplyr::filter(`TGFb-high-isDiffPeak` == TRUE) %>% 
    dplyr::select(peak)
  
  subsetTib <- subsetTib %>% 
    dplyr::mutate(isRAdiffPeak = peak %in% RA.diffPeaks$peak) %>% 
    dplyr::mutate(isTGFbdiffPeak = peak %in% TGFb.diffPeaks$peak)
  
  #######################
  ##Look at change w/RA##
  
  ##RA high dose
    plotTib <- subsetTib %>%  
    dplyr::select(gene_name, peak, isRAdiffPeak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-high_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  #find median per gene
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()
  
  p1  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#C1292E", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid  - 100kb window") +
    xlab("log2 FC median Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"RAhigh_median.svg")
  ggsave(outputLoc, p1, height = 5, width = 5, units = "in", dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
    R=10000)
  
  boot.ci(boot, type = 'basic')
  
  #find max per gene
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`)
  
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name, peak, isRAdiffPeak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-high_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(max_log2FC_access = max(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% unique()
  
  p2  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#C1292E", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Max log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid  - 100 Kb window") +
    xlab("log2 FC max Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"RAhigh_max.svg")
  ggsave(outputLoc, p2, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "max_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  ##RA med dose
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name, peak, isRAdiffPeak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-med-avgNormFragmentCounts`, `RA-med_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-med_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  #find median per gene
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(median_log2FC_access = median(log2(`RA-med-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`RA-med_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()
  
  p3  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#DF686C", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median log2FC Acessibility vs. log2FC Expression in med dose Retinoic acid") +
    xlab("log2 FC median Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"RAmed_median.svg")
  ggsave(outputLoc, p3, height = 5, width = 5, units = "in", dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  #find max per gene
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name, peak, isRAdiffPeak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-med-avgNormFragmentCounts`, `RA-med_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-med_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(max_log2FC_access = max(log2(`RA-med-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`RA-med_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% unique()
  
  p4  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#DF686C", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Max log2FC Acessibility vs. log2FC Expression in med dose Retinoic acid") +
    xlab("log2 FC max Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"RAmed_max.svg")
  ggsave(outputLoc, p4, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "max_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  ##RA low dose
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name, peak, isRAdiffPeak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-low-avgNormFragmentCounts`, `RA-low_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-low_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  #find median per gene
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(median_log2FC_access = median(log2(`RA-low-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`RA-low_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()
  
  p5  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#EDABAD", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median log2FC Acessibility vs. log2FC Expression in low dose Retinoic acid ") +
    xlab("log2 FC median Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"RAlow_median.svg")
  ggsave(outputLoc, p5, height = 5, width = 5, units = "in", dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  #find max per gene
  
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name, peak, isRAdiffPeak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-low-avgNormFragmentCounts`, `RA-low_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-low_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(max_log2FC_access = max(log2(`RA-low-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`RA-low_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% unique()
  
  p6  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#EDABAD", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Max log2FC Acessibility vs. log2FC Expression in low dose Retinoic acid  ") +
    xlab("log2 FC max Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"RAlow_max.svg")
  ggsave(outputLoc, p6, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "max_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  #######################
  ##Look at change w/TGFb##
  
  ##TGFb high dose
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-high_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  #find median per gene
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(median_log2FC_access = median(log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()
  
  p7  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#091E5F", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose TGFb  - 100 Kb window") +
    xlab("log2 FC median Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"TGFbhigh_median.svg")
  ggsave(outputLoc, p7, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  #find max per gene
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-high_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(max_log2FC_access = max(log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% unique()
  
  p8  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#091E5F", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Max log2FC Acessibility vs. log2FC Expression in high dose TGFb  - 100 Kb window") +
    xlab("log2 FC max Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"TGFbhigh_max.svg")
  ggsave(outputLoc, p8, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "max_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  ##TGFb med dose
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-med-avgNormFragmentCounts`, `TGFb-med_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-med_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  #find median per gene
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(median_log2FC_access = median(log2(`TGFb-med-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`TGFb-med_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()
  
  p9  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#306FD8", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median log2FC Acessibility vs. log2FC Expression in med dose TGFb") +
    xlab("log2 FC median Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"TGFbmed_median.svg")
  ggsave(outputLoc, p9, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  
  #find max per gene
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-med-avgNormFragmentCounts`, `TGFb-med_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-med_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(max_log2FC_access = max(log2(`TGFb-med-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`TGFb-med_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% unique()
  
  p10  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#306FD8", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Max log2FC Acessibility vs. log2FC Expression in med dose TGFb") +
    xlab("log2 FC max Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"TGFbmed_max.svg")
  ggsave(outputLoc, p10, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "max_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  ##TGFb low dose
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-low-avgNormFragmentCounts`, `TGFb-low_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-low_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  
  #find median per gene
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(median_log2FC_access = median(log2(`TGFb-low-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`TGFb-low_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()
  
  p11  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#BCE7FD", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median log2FC Acessibility vs. log2FC Expression in low dose TGFb") +
    xlab("log2 FC median Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"TGFblow_median.svg")
  ggsave(outputLoc, p11, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(median_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')
  
  #find max per gene
  plotTib <- subsetTib %>%  
    dplyr::select(gene_name,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-low-avgNormFragmentCounts`, `TGFb-low_avgTPM`)
  
  #remove peaks from genes where the rna value for both samples is zero
  plotTib <- plotTib %>% 
    dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-low_avgTPM` == 0)
  
  plotTib <- plotTib %>% 
    dplyr::filter(expressionFilter  == FALSE)
  #replace zero values with 1 for calculating fold changes
  plotTib[plotTib == 0] = 1
  
  plotTib <- plotTib %>% 
    group_by(gene_name) %>% 
    dplyr::mutate(max_log2FC_access = max(log2(`TGFb-low-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))
  
  plotTib <- plotTib %>% 
    dplyr::mutate(log2FC_expression = log2(`TGFb-low_avgTPM`/`EtOH-nlDensity_avgTPM`))
  
  plotTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% unique()
  
  p12  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#BCE7FD", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Max log2FC Acessibility vs. log2FC Expression in low dose TGFb") +
    xlab("log2 FC max Accessibility (norm. fragment counts)") +
    ylab("log2 FC Expression ")
  
  outputLoc <- paste0(outputFolder,"TGFblow_max.svg")
  ggsave(outputLoc, p12, width = 5, height = 5, units = 'in', dpi = 300)
  
  #calculate bootstrapped R values
  bootTib <- plotTib %>% 
    dplyr::select(max_log2FC_access, log2FC_expression) %>% 
    na.omit()
  
  boot <- boot(bootTib, statistic = function(bootTib, i){
    cor(bootTib[i, "max_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
  },
  R=10000)
  
  boot.ci(boot, type = 'basic')

}