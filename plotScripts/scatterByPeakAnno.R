theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksFile <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))
  outputFolder     <- here('plots', 'scatter')
} else {
  peaksFile     <- cmdargs[1]
  outputLoc     <- cmdargs[2]
}

#select columns needed
subsetTib <- peaksFile %>% 
  dplyr::select(peak, gene_name, annotation,
                `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,
                `RA-low-avgNormFragmentCounts`, `RA-med-avgNormFragmentCounts`, `RA-high-avgNormFragmentCounts`,
                `TGFb-low-avgNormFragmentCounts`, `TGFb-med-avgNormFragmentCounts`, `TGFb-high-avgNormFragmentCounts`,
                `RA-low_avgTPM`, `RA-med_avgTPM`, `RA-high_avgTPM`,
                `TGFb-low_avgTPM`, `TGFb-med_avgTPM`, `TGFb-high_avgTPM`)

#adjust annotation formatting
subsetTib <- subsetTib %>% 
  tidyr::separate(annotation, into = c("anno", NA), sep = "[(]")

subsetTib <- subsetTib %>% 
  dplyr::mutate(annotation = case_when(anno == "Intron " ~ "gene body",
                                       anno == "Exon " ~ "gene body",
                                       anno == "3' UTR" ~ "gene body",
                                       anno == "5' UTR" ~ "gene body",
                                       anno == "Promoter " ~ "promoter",
                                       anno == "Downstream " ~ "downstream",
                                       TRUE ~ "intergenic"))


#RA high dose
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, annotation, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`)) %>% 
  dplyr::mutate(log2FC_access = log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

#find median per anno  per gene
plotTib <- plotTib %>% 
  group_by(gene_name, annotation) %>% 
  dplyr::mutate(median_log2FC_access = median(log2FC_access))

plotTib <- plotTib %>% 
  dplyr::select(annotation, median_log2FC_access, log2FC_expression) %>% unique()

p1  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", facet.by = "annotation", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = -2.5, label.y = 5)  +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"/RAhigh_peakAnno_median.svg")
ggsave(outputLoc, p1, height = 5, width = 5, units = "in", dpi = 300)

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(annotation == "downstream") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#

bootTib <- plotTib %>% 
  dplyr::filter(annotation == "gene body") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#

bootTib <- plotTib %>% 
  dplyr::filter(annotation == "intergenic") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#

bootTib <- plotTib %>% 
  dplyr::filter(annotation == "promoter") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

####
#TGFb high dose
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, annotation, `EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,`TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`)) %>% 
  dplyr::mutate(log2FC_access = log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

#find median per anno  per gene
plotTib <- plotTib %>% 
  group_by(gene_name, annotation) %>% 
  dplyr::mutate(median_log2FC_access = median(log2FC_access))

plotTib <- plotTib %>% 
  dplyr::select(annotation, median_log2FC_access, log2FC_expression) %>% unique()

p2  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", facet.by = "annotation", alpha = 0.5, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = -2.5, label.y = 5)  +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"/TGFbhigh_peakAnno_median.svg")
ggsave(outputLoc, p2, height = 5, width = 5, units = "in", dpi = 300)

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(annotation == "downstream") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#

bootTib <- plotTib %>% 
  dplyr::filter(annotation == "gene body") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#

bootTib <- plotTib %>% 
  dplyr::filter(annotation == "intergenic") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

#

bootTib <- plotTib %>% 
  dplyr::filter(annotation == "promoter") %>% 
  dplyr::select(median_log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "median_log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

