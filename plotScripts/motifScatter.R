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
  dplyr::select(matches("(ENSEMBL)|(_rna)|(peak)|(_normATAC)|(_avgTPM)|(NFKB)")) 
                

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)


#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`EtOH-nlDensity-avgNormFragmentCounts` = `EtOH-nlDensity-avgNormFragmentCounts` + 1) %>%
  dplyr::mutate(`RA-high-avgNormFragmentCounts` = `RA-high-avgNormFragmentCounts` + 1) %>% 
  
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`RA-high_avgTPM` = `RA-high_avgTPM` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))

plotTib <- plotTib %>% arrange(desc(RARA_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isRARApeak = RARA_numMotifMatches > 0)


p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isRARApeak == TRUE), aes(col = RARA_motifMatchScore), alpha = 0.7, size=0.1) +
  scale_color_gradient(low = "#FCEEEF", high = "#C1292E") +xlim(-3,3) +ylim(-4, 6.5) + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isRARApeak == FALSE)), n = 5000) %>% dplyr::filter(isRARApeak == FALSE), col = "black", alpha = 0.8, size = 0.1) +
  xlim(-3,3) +ylim(-4, 6.5)


p1 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/RARA.svg")
ggsave(plot = p1, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)


#bootstrap R values

#calculate bootstrapped R values
bootTib <-  plotTib %>% 
  dplyr::filter(isRARApeak == TRUE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <-  plotTib %>% 
  dplyr::filter(isRARApeak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')

##HOXA13 motif
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, peak, isRAdiffPeak, `RA-high_isDeGene`,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,
                `RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`, HOXA13_motifMatchScore, HOXA13_numMotifMatches)

# #filter for motif positive
# plotTib <- plotTib %>% 
#   dplyr::filter(RARA_numMotifMatches > 0)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)


#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`EtOH-nlDensity-avgNormFragmentCounts` = `EtOH-nlDensity-avgNormFragmentCounts` + 1) %>%
  dplyr::mutate(`RA-high-avgNormFragmentCounts` = `RA-high-avgNormFragmentCounts` + 1) %>% 
  
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`RA-high_avgTPM` = `RA-high_avgTPM` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))

plotTib <- plotTib %>% arrange(desc(HOXA13_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isHOXA13peak = HOXA13_numMotifMatches > 0)

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isHOXA13peak == TRUE), aes(col = HOXA13_motifMatchScore), alpha = 0.7, size=0.1) +
  scale_color_gradient(low = "#FCEEEF", high = "#C1292E") + ylim(-3, 6.5) + xlim(-3,3.5) + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isHOXA13peak == FALSE)), n = 5000) %>% dplyr::filter(isHOXA13peak == FALSE), col = "black", alpha = 0.8, size = 0.1) +
  ylim(-3, 6.5) + xlim(-3,3.5)


p2 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)
  
outputLoc <- paste0(outputFolder,"/HOXA13.svg")
ggsave(plot = p2, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
 dplyr::filter(isHOXA13peak == TRUE) %>% 
 dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isHOXA13peak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')

##FOXA1 motif
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, peak, isRAdiffPeak, `RA-high_isDeGene`,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,
                `RA-high-avgNormFragmentCounts`, `RA-high_avgTPM`, FOXA1_motifMatchScore, FOXA1_numMotifMatches)

# #filter for motif positive
# plotTib <- plotTib %>% 
#   dplyr::filter(RARA_numMotifMatches > 0)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `RA-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)


#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`EtOH-nlDensity-avgNormFragmentCounts` = `EtOH-nlDensity-avgNormFragmentCounts` + 1) %>%
  dplyr::mutate(`RA-high-avgNormFragmentCounts` = `RA-high-avgNormFragmentCounts` + 1) %>% 
  
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`RA-high_avgTPM` = `RA-high_avgTPM` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))

plotTib <- plotTib %>% arrange(desc(FOXA1_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isFOXA1peak = FOXA1_numMotifMatches > 0)

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isFOXA1peak == TRUE), aes(col = FOXA1_motifMatchScore), alpha = 0.7, size=0.1) +
  scale_color_gradient(low = "#FCEEEF", high = "#C1292E") +ylim(-3,6) +xlim(-3,3) + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isFOXA1peak == FALSE)), n = 5000) %>% dplyr::filter(isFOXA1peak == FALSE), col = "black", alpha = 0.8, size = 0.1) +
 ylim(-3,6) +xlim(-3,3)

p3 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/FOXA1.svg")
ggsave(plot = p3, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(isFOXA1peak == TRUE) %>% 
dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isFOXA1peak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')
#########################
##Look at change w/TGFb##

##TGFb high dose

##SMAD3 motif
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, peak, isTGFbdiffPeak, `TGFb-high_isDeGene`,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,
                `TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`, SMAD3_motifMatchScore, SMAD3_numMotifMatches)

# #filter for motif positive
# plotTib <- plotTib %>% 
#   dplyr::filter(RARA_numMotifMatches > 0)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)


#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`EtOH-nlDensity-avgNormFragmentCounts` = `EtOH-nlDensity-avgNormFragmentCounts` + 1) %>%
  dplyr::mutate(`TGFb-high-avgNormFragmentCounts` = `TGFb-high-avgNormFragmentCounts` + 1) %>% 
  
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`TGFb-high_avgTPM` = `TGFb-high_avgTPM` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))

plotTib <- plotTib %>% arrange(desc(SMAD3_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isSMAD3peak = SMAD3_numMotifMatches > 0)

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isSMAD3peak == TRUE), aes(col = SMAD3_motifMatchScore), alpha = 0.7, size=0.1) +
  scale_color_gradient(low = "#DAE2FC", high = "#091E5F") +ylim(-2.5,3) + xlim(-2,4.5) + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isSMAD3peak == FALSE)), n = 5000) %>% dplyr::filter(isSMAD3peak == FALSE), col = "black", alpha = 0.8, size = 0.1) +
  ylim(-2.5,3) + xlim(-2,4.5)

p4 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/SMAD3.svg")
ggsave(plot = p4, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(isSMAD3peak == TRUE) %>% 
dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isSMAD3peak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')

##SMAD4 motif
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, peak, isTGFbdiffPeak, `TGFb-high_isDeGene`,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,
                `TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`, SMAD4_motifMatchScore, SMAD4_numMotifMatches)

# #filter for motif positive
# plotTib <- plotTib %>% 
#   dplyr::filter(RARA_numMotifMatches > 0)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)


#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`EtOH-nlDensity-avgNormFragmentCounts` = `EtOH-nlDensity-avgNormFragmentCounts` + 1) %>%
  dplyr::mutate(`TGFb-high-avgNormFragmentCounts` = `TGFb-high-avgNormFragmentCounts` + 1) %>% 
  
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`TGFb-high_avgTPM` = `TGFb-high_avgTPM` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))

plotTib <- plotTib %>% arrange(desc(SMAD4_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isSMAD4peak = SMAD4_numMotifMatches > 0)

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isSMAD4peak == TRUE), aes(col = SMAD4_motifMatchScore), alpha = 0.7, size=0.1) +
  scale_color_gradient(low = "#DAE2FC", high = "#091E5F") + ylim(-3,4) + xlim(-2,4) + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isSMAD4peak == FALSE)), n = 5000) %>% dplyr::filter(isSMAD4peak == FALSE), col = "black", alpha = 0.8, size = 0.1) +
  ylim(-3,4) + xlim(-2,4)

p5 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/SMAD4.svg")
ggsave(plot = p5, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(isSMAD4peak == TRUE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isSMAD4peak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')

##SMAD9 motif
plotTib <- subsetTib %>%  
  dplyr::select(gene_name, peak, isTGFbdiffPeak, `TGFb-high_isDeGene`,`EtOH-nlDensity-avgNormFragmentCounts`, `EtOH-nlDensity_avgTPM`,
                `TGFb-high-avgNormFragmentCounts`, `TGFb-high_avgTPM`, SMAD9_motifMatchScore, SMAD9_numMotifMatches)

# #filter for motif positive
# plotTib <- plotTib %>% 
#   dplyr::filter(RARA_numMotifMatches > 0)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = `EtOH-nlDensity_avgTPM` == 0 & `TGFb-high_avgTPM` == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)


#Add 1 normalized fragment count or tpm to avoid dividing by zero
plotTib <- plotTib %>%
  dplyr::mutate(`EtOH-nlDensity-avgNormFragmentCounts` = `EtOH-nlDensity-avgNormFragmentCounts` + 1) %>%
  dplyr::mutate(`TGFb-high-avgNormFragmentCounts` = `TGFb-high-avgNormFragmentCounts` + 1) %>% 
  
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`TGFb-high_avgTPM` = `TGFb-high_avgTPM` + 1)


# plotTib <- plotTib %>%
#   group_by(gene_name) %>%
#   dplyr::mutate(median_log2FC_access = median(log2(`RA-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`)))

plotTib <- plotTib %>%
  dplyr::mutate(log2FC_access = log2(`TGFb-high-avgNormFragmentCounts`/`EtOH-nlDensity-avgNormFragmentCounts`))

plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))

plotTib <- plotTib %>% arrange(desc(SMAD9_motifMatchScore))

plotTib <- plotTib %>% 
  dplyr::mutate(isSMAD9peak = SMAD9_numMotifMatches > 0)

p.motif <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = plotTib %>% dplyr::filter(isSMAD9peak == TRUE), aes(col = SMAD9_motifMatchScore), alpha = 0.7, size=0.1) +
  scale_color_gradient(low = "#DAE2FC", high = "#091E5F") + ylim(-2.5,4) + xlim(-2,4) + theme(legend.position = "")

p.not <- ggplot(plotTib, aes(x = log2FC_access, log2FC_expression)) +
  geom_point(data = slice_sample((plotTib %>% dplyr::filter(isSMAD9peak == FALSE)), n = 5000) %>% dplyr::filter(isSMAD9peak == FALSE), col = "black", alpha = 0.8, size = 0.1) +
  ylim(-2.5,4) + xlim(-2,4)

p6 <- ggarrange(p.motif,p.not, ncol = 2, nrow = 1)

outputLoc <- paste0(outputFolder,"/SMAD9.svg")
ggsave(plot = p6, filename = outputLoc, width = 16, height = 10, units = "in", dpi = 300)

#bootstrap R values

#calculate bootstrapped R values
bootTib <- plotTib %>% 
  dplyr::filter(isSMAD9peak == TRUE) %>% 
dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=10000)

boot.ci(boot, type = 'basic')

bootTib <- plotTib %>% 
  dplyr::filter(isSMAD9peak == FALSE) %>% 
  dplyr::select(log2FC_access, log2FC_expression) %>% 
  na.omit()

boot <- boot(bootTib, statistic = function(bootTib, i){
  cor(bootTib[i, "log2FC_access"], bootTib[i, "log2FC_expression"], method = "pearson")
},
R=1000)

boot.ci(boot, type = 'basic')
