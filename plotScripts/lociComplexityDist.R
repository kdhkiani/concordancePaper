theme_set(theme_classic())

##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib   <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))
  outputLoc <- here('plots')
} else {
  bigtib    <- read_tsv(cmdargs[1])
  outputLoc  <- cmdargs[3]
}


##massage and prepare data from bigtib

#pull out the relevant values 
subsetTib <- bigtib %>% 
  dplyr::select(gene_name, peak,
                `TGFb-high_avgTPM`, `RA-high_avgTPM`, `EtOH-nlDensity-avgNormFragmentCounts`, `RA-high-avgNormFragmentCounts`,
                `RA-med-avgNormFragmentCounts`, `RA-low-avgNormFragmentCounts`, `TGFb-high-avgNormFragmentCounts`,
                `TGFb-med-avgNormFragmentCounts`, `TGFb-low-avgNormFragmentCounts`) %>% 
  tidyr::separate(peak, into = c(NA, "startLocs", "endLocs"), sep = ":") %>% 
  dplyr::mutate(peakWidth = abs(as.numeric(endLocs) - as.numeric(startLocs)) )


#annotate max peak value between RA/etOH & TGFb/etOH and if greater than or less than 30
subsetTib <- subsetTib %>% 
  rowwise() %>% 
  dplyr::mutate(RA.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`RA-high-avgNormFragmentCounts`,
                                `RA-med-avgNormFragmentCounts`, `RA-low-avgNormFragmentCounts`) > 30) %>% 
  dplyr::mutate(TGFb.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`TGFb-high-avgNormFragmentCounts`,
                                  `TGFb-med-avgNormFragmentCounts`, `TGFb-low-avgNormFragmentCounts`) > 30)
#subset on peaks that pass filter
RA.peaks <- subsetTib %>% dplyr::filter(RA.maxVal == TRUE)
TGFb.peaks <- subsetTib %>% dplyr::filter(TGFb.maxVal == TRUE)

#get counts for number of peaks per gene
RA.peakCounts <- RA.peaks %>% dplyr::count(gene_name)
TGFb.peakCounts <- TGFb.peaks %>% dplyr::count(gene_name)

#combine peak distributions and make plot
plotTib <- RA.peakCounts %>% 
  dplyr::select(n) %>% 
  dplyr::rename(value = n) %>% 
  dplyr::mutate(signal = "RA")

tmp <- TGFb.peakCounts %>% 
  dplyr::select(n) %>% 
  dplyr::rename(value = n) %>% 
  dplyr::mutate(signal = "TGFb")

plotTib <- bind_rows(plotTib, tmp)

RA.tertile <- quantile(RA.peakCounts$n, probs = seq(0, 1,by=1/3))[3]
TGFb.tertile <- quantile(TGFb.peakCounts$n, probs = seq(0, 1,by=1/3))[3]

RA.tertile == TGFb.tertile

p1 <- ggdensity(plotTib, x = "value",
                  add = "median", color = "signal", fill = "signal",
                  palette = c("#C1292E", "#091E5F")) + xlim(0, 20) + geom_vline(xintercept = RA.tertile)

outputFile <- paste0(outputLoc, "/complexityDensity.svg")
ggsave(filename = outputFile, plot = p1)


#annotate gene complexity type based on defintions from gonzalez et al, 2015
subsetTib <- left_join(subsetTib, RA.peakCounts) %>% 
  dplyr::rename(RA.nPeaks = n) %>% 
  left_join(TGFb.peakCounts) %>% 
  dplyr::rename(TGFb.nPeaks = n) %>% 
  dplyr::mutate(RA_lociComplexity = case_when(RA.nPeaks > 7 ~ "high",
                                           TRUE ~ "low")) %>% 
  dplyr::mutate(TGFb_lociComplexity = case_when(TGFb.nPeaks > 7 ~ "high",
                                              TRUE ~ "low"))

#pull data to plot with 
plotTib <- subsetTib %>% 
  dplyr::select(gene_name, `RA-high_avgTPM`, `TGFb-high_avgTPM`, RA_lociComplexity, TGFb_lociComplexity) %>% 
  unique() %>% 
  dplyr::mutate(RA_log2tpm = log2(`RA-high_avgTPM`)) %>% 
  dplyr::mutate(TGFb_log2tpm = log2(`TGFb-high_avgTPM`)) %>% 
  dplyr::select(-`RA-high_avgTPM`, -`TGFb-high_avgTPM`)

#make it long format after much struggle for plotting
plotTib <- plotTib %>% tidyr::gather(key = "signal", value = "log2tpm", -gene_name, -RA_lociComplexity, -TGFb_lociComplexity) %>% 
  tidyr::gather(key = "signal2", value = "lociComplexity", -gene_name, -signal, -log2tpm) %>% 
  dplyr::select(-signal2) %>% 
  tidyr::separate(signal, c("signal", NA), sep = "_") %>% 
  dplyr::mutate(lociComplexity = factor(lociComplexity, levels = c("low", "high")))

p2 <-  plotTib %>% 
  ggplot(aes(x=lociComplexity, y=log2tpm, color = lociComplexity, fill = lociComplexity)) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width=0, justification = -0.3, point_colour=NA) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(size = 0.2, alpha = 0.4, position = position_jitter(seed=1, width = 0.1)) +
  scale_color_manual(values = c("#1B9aaa","#56e39f")) + 
  scale_fill_manual(values = c("#1B9aaa","#56e39f")) + 
  facet_wrap(~signal) +theme(legend.position = "")

ggsave(paste0(outputLoc, "/complexityExpression.svg") , plot = p2)

###now lets look at peakwidths
plotTib <- subsetTib %>% 
  dplyr::select(peakWidth, RA_lociComplexity)

p3 <- plotTib %>% 
  dplyr::mutate(RA_lociComplexity= factor(RA_lociComplexity, levels = c("low", "high"))) %>% 
  ggplot(aes(x=peakWidth, group=RA_lociComplexity, fill=RA_lociComplexity)) +
  geom_density(adjust=1.5) +
  facet_wrap(~ RA_lociComplexity) +
  scale_x_continuous(limits =c(145,500), breaks = seq(150, 500, by = 50)) +
  scale_fill_manual(values = c("#1B9aaa","#56e39f")) + theme(legend.position = "")

ggsave(paste0(outputLoc, "/complexity/complexityPeak.svg"), plot = p3, width = 5, height = 5, units = "in", dpi = 300)


