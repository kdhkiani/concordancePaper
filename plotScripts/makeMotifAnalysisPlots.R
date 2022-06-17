register(MulticoreParam(4, progressbar = TRUE))
theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCountsDiffPeaks <- readRDS(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.rds'))
  selected.PWM.objects    <- readRDS(here('extractedData', 'most_variable_motifs_Robject.rds'))
  outputLoc               <- here('plots', 'validation')
} else {
  fragmentCountsDiffPeaks <- readRDS(cmdargs[1])
  selected.PWM.objects    <- readRDS(cmdargs[2])
  outputLoc       <- cmdargs[3]
}


tf.factor.order <- c("RARA", "SMAD3", "SMAD4", "SMAD9", 
                     "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", 
                     "BACH1",  "BACH2", "BATF", "FOXA1", "FOXA2", "FOXA3", "FOXC2", "FOXD3", 
                     "HOXA13", "HOXB13", "HOXC10", "HOXC12", "HOXC13", "HOXD13", 
                     "NFKB1", "REL", "RELA", "SMARCC1", "CDX1",  "CDX2", 
                     "CTCF", "NFE2", "NFE2L2", "MAFF", "MAFK", "BCL11A", "BCL11B", 
                     "GRHL1", "SPI1", "SPIB", "SPIC", 
                     "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5", "ELK4", "ETS2")

n.bootstrap.samples <- 1000

#use the selected PWM objects to test for motif deviations in a subset of the data: controls and TGFB only, controls and RA only
makeDeviationScorePlotTibs <- function(fragmentCountsDiffPeaks, motif_ix) {
  allSampleNames <- colnames(fragmentCountsDiffPeaks)
  etohSampleNames <- allSampleNames[grepl("EtOH", allSampleNames)]
  raSampleNames <- allSampleNames[grepl("-RA-", allSampleNames) & !grepl("-TGFb-", allSampleNames)]
  tgfbSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & !grepl("-RA-", allSampleNames)]
  
  # calculate deviation scores looking just the RA condition and how it compares to the EtOH condition
  fragCountsRAandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, raSampleNames)]
  devRA <- computeDeviations(object = fragCountsRAandControls, annotations = motif_ix)
  variabilityRA <- computeVariability(devRA)
  tvarRA <- as.tibble(variabilityRA)
  devScoresEtOHconds <- (assays(devRA)$deviations[, 1:9])
  devScoresRAconds   <- (assays(devRA)$deviations[, 10:18])
  avgDevEtOH <- rowMeans(devScoresEtOHconds)
  avgDevRA   <- rowMeans(devScoresRAconds)
  devScoreRA <- (avgDevRA + 1) / (avgDevEtOH + 1) - 1
  tfnames <- names(avgDevEtOH)
  tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
  n_tfs = length(tfnames)
  raDevTib <- tibble(tf_name = factor(tfnames, levels = tf.factor.order), dev_score = devScoreRA, cond = rep("RA", n_tfs))
  
  # calculate deviation scores looking just the TGFb condition and how it compares to the EtOH condition
  fragCountsTGFbandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, tgfbSampleNames)]
  devTGFb <- computeDeviations(object = fragCountsTGFbandControls, annotations = motif_ix)
  variabilityTGFb <- computeVariability(devTGFb)
  tvarTGFb <- as.tibble(variabilityTGFb)
  devScoresEtOHconds <- (assays(devTGFb)$deviations[, 1:9])
  devScoresTGFbconds   <- (assays(devTGFb)$deviations[, 10:18])
  avgDevEtOH <- rowMeans(devScoresEtOHconds)
  avgDevTGFb   <- rowMeans(devScoresTGFbconds)
  devScoreTGFb <- (avgDevTGFb + 1) / (avgDevEtOH + 1) - 1
  tfnames <- names(avgDevEtOH)
  tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
  n_tfs = length(tfnames)
  TGFbDevTib <- tibble(tf_name = factor(tfnames, levels = tf.factor.order), dev_score = devScoreTGFb, cond = rep("TGFb", n_tfs))
  
  combTib <- rbind(raDevTib, TGFbDevTib) 
  
  return(combTib)  
}

fragmentCountsDiffPeaks <- addGCBias(fragmentCountsDiffPeaks, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ix <- matchMotifs(selected.PWM.objects, fragmentCountsDiffPeaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devScoresTib <- makeDeviationScorePlotTibs(fragmentCountsDiffPeaks, motif_ix)

# add bootstrap confidence intervals
set.seed(0)
n.peaks <- nrow(fragmentCountsDiffPeaks)
bootstrap.dev.scores.tib <- NULL
for (ii in 1:n.bootstrap.samples) {
  this.sample.inds <- sample(1:n.peaks, n.peaks, replace = TRUE)
  fragCtsBootstrapSample <- fragmentCountsDiffPeaks[this.sample.inds, ]
  this.motif_ix <- motif_ix[this.sample.inds, ]
  this.devScoreTib <- makeDeviationScorePlotTibs(fragCtsBootstrapSample, this.motif_ix)
  this.devScores <- this.devScoreTib$dev_score
  bootstrap.dev.scores.tib <- rbind(bootstrap.dev.scores.tib, this.devScores)
  
  print(sprintf("%d of %d bootstrap samples complete for motif enrichment by TF plot", ii, n.bootstrap.samples))
}
lower.bootstrap.quantiles <- c()
upper.bootstrap.quantiles <- c()
for (ii in 1:ncol(bootstrap.dev.scores.tib)) {
  lower.bootstrap.quantiles <- c(lower.bootstrap.quantiles, quantile(bootstrap.dev.scores.tib[,ii], .05))
  upper.bootstrap.quantiles <- c(upper.bootstrap.quantiles, quantile(bootstrap.dev.scores.tib[,ii], .95))
}
devScoresTib[["bootstrap_ci_lower"]] <- 2 * devScoresTib$dev_score - upper.bootstrap.quantiles
devScoresTib[["bootstrap_ci_upper"]] <- 2 * devScoresTib$dev_score - lower.bootstrap.quantiles


# devscore.plot.list <- list()
# colorscheme <- c('dark green', 'blue', 'dark orange')
# for (condname in c("RA", "TGFb", "Both")) {
#   p <- devScoresTib %>%
#     filter(cond == condname) %>%
#     ggplot(aes(tf_name, dev_score, ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper)) + 
#     geom_bar(stat="identity", fill = colorscheme[length(devscore.plot.list) + 1]) + 
#     geom_errorbar(width = 0) +
#     geom_hline(yintercept = 0) +
#     ylab(paste0("dev score ", condname)) + 
#     xlab("") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5)) 
#   devscore.plot.list[[length(devscore.plot.list) + 1]] <- p
# }

plotTib <- devScoresTib %>% 
  tidyr::drop_na() %>%
  dplyr::filter(tf_name %in% c("RARA", "SMAD3", "SMAD4", "SMAD9", 
                               "HOXA13", "HOXB13", "HOXC10", "HOXC12", "HOXC13", "HOXD13", 
                               "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2"))

p1<- plotTib %>% 
  dplyr::mutate(tf_name = factor(tf_name, levels = c("RARA", "SMAD3", "SMAD4", "SMAD9", 
                                                     "HOXA13", "HOXB13", "HOXC10", "HOXC12", "HOXC13", "HOXD13", 
                                                     "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2"))) %>% 
  ggplot(aes(x=tf_name, y=dev_score, ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper, fill = cond)) +
  geom_bar(stat="identity") + 
  geom_errorbar(width = 0) + theme(axis.text.x = element_text(angle = 45)) +
  facet_grid(cond ~ .) +
  scale_fill_manual(values = c("#C1292E", "#091E5F"))

ggsave(filename = paste0(outputLoc,"/motifDevScore.svg"), plot = p1)
  
  
  
