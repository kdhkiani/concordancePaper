theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigTib        <- read_tsv(here('extractedData', 'Jurida_Weiterer', 'jw_genePeak_1to1mapping_mergedDist50.tsv'))#karun's input
  peakAnno      <- read_tsv(here('extractedData', 'Jurida_Weiterer','ATACseq', 'jw_AtacPeaksMergedSummitWindows_mergedist50.annotated.tsv'))#karun's input
  outputFolder  <- here("plots", "scatter")
} else {
  bigTib        <- read_tsv(cmdargs[1])
  peakAnno      <- read_tsv(cmdargs[2])
  outputFolder   <- cmdargs[3]
}

set.seed(2022)

peakAnno <- peakAnno %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":"))

subsetTib <- dplyr::select(bigTib, matches("(ENSEMBL)|(_rna)|(peak)|(_normATAC)|(_avgTPM)")) 

subsetTib <- left_join(subsetTib, peakAnno)

#######################
##Look at change w/IL1alpha##

##NFKB1 motif
plotTib <- subsetTib %>%  
  dplyr::select(matches("(ENSEMBL)|(_rna)|(peak)|(_normATAC)|(_avgTPM)|(NFKB1)")) 
                

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `untreated_rna` == 0 & `IL1a_1hr_rna` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)

#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`untreated_normATAC` = `untreated_normATAC` + 1) %>%
  dplyr::mutate(`IL1a_1hr_normATAC` = `IL1a_1hr_normATAC` + 1) %>% 
  
  dplyr::mutate(`untreated_rna` = `untreated_rna` + 1) %>% 
  dplyr::mutate(`IL1a_1hr_rna` = `IL1a_1hr_rna` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`IL1a_1hr_normATAC`/`untreated_normATAC`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`IL1a_1hr_rna`/`untreated_rna`))

plotTib <- plotTib %>% arrange(desc(NFKB1_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isNFKB1peak = NFKB1_motifMatchScore > 0)


p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isNFKB1peak == TRUE), aes(col = NFKB1_motifMatchScore), alpha = 0.7, size=0.5) +
  scale_color_gradient(low = "#FEF9EC", high = "#EFC22E") + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isNFKB1peak == FALSE)), n = 5000) %>% dplyr::filter(isNFKB1peak == FALSE), col = "black", alpha = 0.8, size = 0.5) 


p1 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/NFKB1.svg")
ggsave(plot = p1, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)


#bootstrap R values

#calculate bootstrapped R values
bootTib <-  plotTib %>% 
  dplyr::filter(isNFKB1peak == TRUE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <-  plotTib %>% 
  dplyr::filter(isNFKB1peak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')

##REL motif
plotTib <- subsetTib %>%  
  dplyr::select(matches("(ENSEMBL)|(_rna)|(peak)|(_normATAC)|(_avgTPM)|(REL)")) 

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `untreated_rna` == 0 & `IL1a_1hr_rna` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)

#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`untreated_normATAC` = `untreated_normATAC` + 1) %>%
  dplyr::mutate(`IL1a_1hr_normATAC` = `IL1a_1hr_normATAC` + 1) %>% 
  
  dplyr::mutate(`untreated_rna` = `untreated_rna` + 1) %>% 
  dplyr::mutate(`IL1a_1hr_rna` = `IL1a_1hr_rna` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`IL1a_1hr_normATAC`/`untreated_normATAC`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`IL1a_1hr_rna`/`untreated_rna`))

plotTib <- plotTib %>% arrange(desc(REL_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isRELpeak = REL_motifMatchScore > 0)

#plot

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isRELpeak == TRUE), aes(col = REL_motifMatchScore), alpha = 0.7, size=0.5) +
  scale_color_gradient(low = "#FEF9EC", high = "#EFC22E") + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isRELpeak == FALSE)), n = 5000) %>% dplyr::filter(isRELpeak == FALSE), col = "black", alpha = 0.8, size = 0.5) 


p2 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)
  
outputLoc <- paste0(outputFolder,"/REL.svg")
ggsave(plot = p2, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
 dplyr::filter(isRELpeak == TRUE) %>% 
 dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isRELpeak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')


plotTib <- plotTib %>% arrange(desc(RELA_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isRELApeak = RELA_motifMatchScore > 0)

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isRELApeak == TRUE), aes(col = RELA_motifMatchScore), alpha = 0.7, size=0.5) +
  scale_color_gradient(low = "#FEF9EC", high = "#EFC22E") + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isRELApeak == FALSE)), n = 5000) %>% dplyr::filter(isRELApeak == FALSE), col = "black", alpha = 0.8, size = 0.5) 

p3 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/RELA.svg")
ggsave(plot = p3, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(isRELApeak == TRUE) %>% 
dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isRELApeak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')
