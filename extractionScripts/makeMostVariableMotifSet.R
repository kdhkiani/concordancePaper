register(MulticoreParam(4, progressbar = TRUE))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCounts <- read_rds(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.rds'))
  outputFile <- here("extractedData", "most_variable_motifs_Robject.rds")
} else {
  fragmentCounts <- read_rds(cmdargs[1])
  outputFile <- cmdargs[2]
}

# select the most variable motifs from the data set (selected by finding "natural looking cutoff" in the full set of curated motifs)
fragmentCounts <- addGCBias(fragmentCounts, 
                            genome = BSgenome.Hsapiens.UCSC.hg38)

set.seed(2021)
data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper
motif_ix <- matchMotifs(motifSet, fragmentCounts,
                        genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object = fragmentCounts, annotations = motif_ix)
variability <- computeVariability(dev)
p <- plotVariability(variability, use_plotly = FALSE)
top_n_motifs_to_keep <- 150  # this was chosen by manually looking at the plot below to see where the natural break was
p <- p + geom_vline(xintercept = top_n_motifs_to_keep)
tvar <- as_tibble(variability)

TF.names.variabilitycutoff <- arrange(tvar, -variability)$name[1:top_n_motifs_to_keep]
TF.names.manual.additions  <- tvar %>% filter(name %in% manually_include_these_motifs) %>% pull("name")
TF.names.variabilitycutoff.plus.manual.additions <- union(TF.names.variabilitycutoff, TF.names.manual.additions)
PWM.object.indices.bool <- sapply(strsplit(names(motifSet), "_"), function (x) x[3]) %in% TF.names.variabilitycutoff.plus.manual.additions
PWM.object.indices.num  <- which(PWM.object.indices.bool)
selected.PWM.objects <- motifSet[PWM.object.indices.num]
saveRDS(selected.PWM.objects, file = outputFile)
################################################################################################################  