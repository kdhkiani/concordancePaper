cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib           <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))
  intron.bigtib    <- read_tsv(here('extractedData', 'intron_genePeak_1to1mapping_mergedDist50.tsv'))
  outputFolder     <- here('plots', 'scatter', '')
} else {
  bigtib          <- cmdargs[1]
  intron.bigtib   <- cmdargs[2]
  outputFolder    <- cmdargs[3]

}

set.seed(2022)

####Reformat data to get raw RNA counts and quantile normalize them###

##start with intron data
intron.bigtib <- intron.bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(EtOH_nlDensity_intron_avgCounts = mean(c(`05-EtOH-nlDensity_intron`, `15-EtOH-nlDensity_intron`, `27-EtOH-nlDensity_intron`))) %>%
  dplyr::mutate(RA_high_intron_avgCounts = mean(c(`06-RA-high_intron`, `23-RA-high_intron`, `32-RA-high_intron`))) %>% 
  dplyr::mutate(TGFb_high_intron_avgCounts = mean(c(`10-TGFb-high_intron`, `22-TGFb-high_intron`, `31-TGFb-high_intron`)))
  
intron.subsetTib <- intron.bigtib %>% 
  dplyr::select(ensg, peak, ends_with("intron_avgCounts"), ends_with("avgNormFragmentCounts")) %>% 
  dplyr::select(-contains("-and-"), -contains("-med-"), -contains("-low-"), -contains("-highDensity-"), -contains("-halfDensity-"))
  
intron_mat <- as.matrix(intron.subsetTib[,c(3:5)])
quant.norm.intron.mat <- normalize.quantiles(intron_mat)
colnames(quant.norm.intron.mat) <- paste0(colnames(intron_mat), "_quantNorm")
quant.norm.intron.tib <- as_tibble(quant.norm.intron.mat)

intron.subsetTib <- intron.subsetTib %>% 
  dplyr::bind_cols(quant.norm.intron.tib)

##now with 'normal' exon-based read counts
subsetTib <- bigtib %>% 
  dplyr::select(ensg, peak, ends_with("_avgCounts"), ends_with("avgNormFragmentCounts")) %>% 
  dplyr::select(-contains("-and-"), -contains("-med_"), -contains("-low_"), -contains("-highDensity_"), -contains("-halfDensity_"))

exon_mat <- as.matrix(subsetTib[,c(3:5)])
quant.norm.exon.mat <- normalize.quantiles(exon_mat)
colnames(quant.norm.exon.mat) <- paste0(colnames(exon_mat), "_quantNorm")
quant.norm.exon.tib <- as_tibble(quant.norm.exon.mat)

subsetTib <- subsetTib %>% 
  dplyr::bind_cols(quant.norm.exon.tib)

#######################
##Look at change w/RA##

##RA high dose - exon
plotTib <- subsetTib %>%  
  dplyr::select(ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgCounts_quantNorm` ,`RA-high-avgNormFragmentCounts`, `RA-high_avgCounts_quantNorm`)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgCounts_quantNorm` == 0 & `RA-high_avgCounts_quantNorm` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


#find median per gene
plotTib <- plotTib %>% 
  group_by(ensg) %>% 
  dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`RA-high_avgCounts_quantNorm`/`EtOH-nlDensity_avgCounts_quantNorm`))

plotTib <- plotTib %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()

p1  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#C1292E", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid  - Exon counts") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"exonCounts_RAhigh_median.svg")
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

##RA high dose - intron
plotTib <- intron.subsetTib %>%  
  dplyr::select(ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH_nlDensity_intron_avgCounts_quantNorm` ,`RA-high-avgNormFragmentCounts`, `RA_high_intron_avgCounts_quantNorm`)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH_nlDensity_intron_avgCounts_quantNorm` == 0 & `RA_high_intron_avgCounts_quantNorm` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


#find median per gene
plotTib <- plotTib %>% 
  group_by(ensg) %>% 
  dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`RA_high_intron_avgCounts_quantNorm`/`EtOH_nlDensity_intron_avgCounts_quantNorm`))

plotTib <- plotTib %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()

p2  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#C1292E", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid - intron counts") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"intronCounts_RAhigh_median.svg")
ggsave(outputLoc, p2, height = 5, width = 5, units = "in", dpi = 300)

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#########################
##Look at change w/TGFb##

##TGFb high dose - exon
plotTib <- subsetTib %>%  
  dplyr::select(ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgCounts_quantNorm` ,`TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgCounts_quantNorm`)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgCounts_quantNorm` == 0 & `TGFb-high_avgCounts_quantNorm` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


#find median per gene
plotTib <- plotTib %>% 
  group_by(ensg) %>% 
  dplyr::mutate(median_log2FC_access = median(log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgCounts_quantNorm`/`EtOH-nlDensity_avgCounts_quantNorm`))

plotTib <- plotTib %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()

p3  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#091E5F", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid  - Exon counts") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"exonCounts_TGFbhigh_median.svg")
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

##TGFb high dose - intron
plotTib <- intron.subsetTib %>%  
  dplyr::select(ensg, peak, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH_nlDensity_intron_avgCounts_quantNorm` ,`TGFb-high-avgNormFragmentCounts`, `TGFb_high_intron_avgCounts_quantNorm`)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH_nlDensity_intron_avgCounts_quantNorm` == 0 & `TGFb_high_intron_avgCounts_quantNorm` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


#find median per gene
plotTib <- plotTib %>% 
  group_by(ensg) %>% 
  dplyr::mutate(median_log2FC_access = median(log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`TGFb_high_intron_avgCounts_quantNorm`/`EtOH_nlDensity_intron_avgCounts_quantNorm`))

plotTib <- plotTib %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% unique()

p4  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", color = "#091E5F", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
  ggtitle("Median log2FC Acessibility vs. log2FC Expression in high dose Retinoic acid - intron counts") +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"intronCounts_TGFbhigh_median.svg")
ggsave(outputLoc, p4, height = 5, width = 5, units = "in", dpi = 300)

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')


