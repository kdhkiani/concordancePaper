cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigTib     <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))
  diffPeaks  <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.tsv'))
} else {
  bigTib  <- read_tsv(cmdargs[1])
  diffPeaks <- read_tsv(cmdargs[2])
}

subsetTib <- bigTib %>% 
  dplyr::select(ensg, gene_name, peak, `RA-high_isDeGene`, `TGFb-high_isDeGene`)

diffPeaks <- diffPeaks %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":"))

RA.diffPeaks <- diffPeaks %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE)

TGFb.diffPeaks <- diffPeaks %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE)

subsetTib <- subsetTib %>% 
  dplyr::mutate(isRAdiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::mutate(isTGFbdiffPeak = peak %in% TGFb.diffPeaks$peak)

total.genes <- length(subsetTib$gene_name %>% unique())

###

RA.genes <- subsetTib %>% 
  dplyr::filter(`RA-high_isDeGene` == 1)

RA.genes <- unique(RA.genes$gene_name)
numRA.genes <- length(RA.genes)

RA.peaks <- subsetTib %>% 
  dplyr::filter(isRAdiffPeak == TRUE)

RA.peaks <- unique(RA.peaks$gene_name)
numRA.peaks <- length(RA.peaks)

numRAoverlap <- length(intersect(RA.genes, RA.peaks))

nonRA.genes <- subsetTib %>% 
  dplyr::filter(`RA-high_isDeGene` == 0)
nonRA.genes <- unique(nonRA.genes$gene_name)

nonRA.peaks <- subsetTib %>% 
  dplyr::filter(isRAdiffPeak == FALSE)
nonRA.peaks <- unique(nonRA.peaks$gene_name)

numNonOverlap <- length(intersect(nonRA.genes, nonRA.peaks)) 

#make contingency table for fisher's exact test
A <- numRAoverlap
B <- numRA.genes - numRAoverlap
C <- numRA.peaks - numRAoverlap
D <- numNonOverlap

contingencyMatrix <- matrix(c(A,C, B,D), nrow = 2, ncol = 2)

fisher.test(contingencyMatrix, alternative = "greater")


###

TGFb.genes <- subsetTib %>% 
  dplyr::filter(`TGFb-high_isDeGene` == 1)

TGFb.genes <- unique(TGFb.genes$gene_name)
numTGFb.genes <- length(TGFb.genes)

TGFb.peaks <- subsetTib %>% 
  dplyr::filter(isTGFbdiffPeak == TRUE)

TGFb.peaks <- unique(TGFb.peaks$gene_name)
numTGFb.peaks <- length(TGFb.peaks)

numTGFboverlap <- length(intersect(TGFb.genes, TGFb.peaks))

nonTGFb.genes <- subsetTib %>% 
  dplyr::filter(`TGFb-high_isDeGene` == 0)
nonTGFb.genes <- unique(nonTGFb.genes$gene_name)

nonTGFb.peaks <- subsetTib %>% 
  dplyr::filter(isTGFbdiffPeak == FALSE)
nonTGFb.peaks <- unique(nonTGFb.peaks$gene_name)

numNonOverlap <- length(intersect(nonTGFb.genes, nonTGFb.peaks)) 

#make contingency table for fisher's exact test
A <- numTGFboverlap
B <- numTGFb.genes - numTGFboverlap
C <- numTGFb.peaks - numTGFboverlap
D <- numNonOverlap

contingencyMatrix <- matrix(c(A,C, B,D), nrow = 2, ncol = 2)

fisher.test(contingencyMatrix, alternative = "greater")
