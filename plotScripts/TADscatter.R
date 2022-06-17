deseqout <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
atacFragCounts <- read_tsv(here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.annotated.tsv'))

tss_TADtib <- read_tsv(here('extractedData', 'tss_TADkey.tsv'), col_names = FALSE)
peak_TADtib <- read_tsv(here('extractedData', 'allPeaks_TADkey.tsv'), col_names = FALSE)
tss <- read_tsv(here('refs', 'EnsgToHg38CanonicalTssMapping.tsv'))
outputFolder <- here('plots', 'scatter', '')

set.seed(2022)

deseqout <- deseqout %>% 
  dplyr::select(matches("(gene_name)|(ensg)|(_avgTPM)")) %>% 
  dplyr::select(-contains("highDensity")) %>%
  dplyr::select(-contains("halfDensity")) %>% 
  dplyr::select(-contains("-and-"))

atacFragCounts <- atacFragCounts %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":")) %>% 
  dplyr::select(matches("(peak)|(avgNormFragmentCounts)")) %>% 
  dplyr::select(-contains("highDensity")) %>%
  dplyr::select(-contains("halfDensity")) %>% 
  dplyr::select(-contains("-and-"))
    
tss_TADtib <- tss_TADtib %>% 
  dplyr::mutate(transcript = paste(X1, X2, X3, sep=":")) %>% 
  dplyr::mutate(TAD = paste(X4, X5, X6, sep=":")) %>%
  dplyr::select(transcript, TAD)

tss <- tss %>% 
  dplyr::mutate(transcript = paste(chrom, tx_start, tx_end, sep=":")) %>% 
  dplyr::rename(ensg = gene_id)

transcript.key <- left_join(tss, tss_TADtib) %>% 
  drop_na()

gene <- left_join(deseqout, transcript.key) %>% 
  drop_na()

peak_TADtib <- peak_TADtib %>% 
  dplyr::mutate(peak = paste(X1, X2, X3, sep=":")) %>% 
  dplyr::mutate(TAD = paste(X4, X5, X6, sep=":")) %>%
  dplyr::select(peak, TAD)

peak <- left_join(atacFragCounts, peak_TADtib) %>% 
  drop_na()

bigTib <- left_join(gene, peak) %>% 
  drop_na()


#######################
##Look at change w/RA##

##RA high dose
plotTib <- bigTib %>%  
  dplyr::select(gene_name, ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`, TAD)

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
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid  - all peaks within gene TAD") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"RAhigh_median_TAD.svg")
ggsave(outputLoc, p1, width = 5, height = 5, units = 'in', dpi = 300)

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
plotTib <- bigTib %>%  
  dplyr::select(gene_name, ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`, TAD)

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
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid  - all peaks within gene TAD") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"RAhigh_max_TAD.svg")
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

##########################
##Look at change w/TGFb##

##TGFb high dose
plotTib <- bigTib %>%  
  dplyr::select(gene_name, ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`, TAD)

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

p3  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#091E5F", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = 0, label.y = 5) + stat_regline_equation(label.x = 0,label.y = 6) + 
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose TGFb  - All Peaks within Gene TAD") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression")

outputLoc <- paste0(outputFolder,"TGFbhigh_median_TAD.svg")
ggsave(outputLoc, p3, width = 5, height = 5, units = 'in', dpi = 300)

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
plotTib <- bigTib %>%  
  dplyr::select(gene_name, ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`, TAD)

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

p4  <- ggscatter(plotTib, x = "max_log2FC_access", y = "log2FC_expression", color = "#091E5F", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = 0, label.y = 6) + stat_regline_equation(label.x = 0,label.y = 7) + 
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose TGFb  - All Peaks within Gene TAD") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression")

outputLoc <- paste0(outputFolder,"TGFbhigh_max_TAD.svg")
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
