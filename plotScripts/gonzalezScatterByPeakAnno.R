theme_set(theme_classic())

##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksanno   <- read_tsv(here('extractedData', 'Gonzalez2015', 'peaksTable.tsv'))
  peaksdat    <- read_tsv(here('extractedData', 'Gonzalez2015', 'DNaseCnts.tsv'))
  bedFile     <- here('extractedData', 'Gonzalez2015', 'gonzalezPeaks.bed')
  deseqtib    <- read_tsv(here('extractedData', 'Gonzalez2015', 'RNAseqCnts.tsv'))
  outputFolder  <- here('plots','scatter')
} else {
  peaksanno   <- read_tsv(cmdargs[1])
  peaksdat    <- read_tsv(cmdargs[2])
  bedFile     <- read_tsv(cmdargs[3])
  deseqtib    <- read_tsv(cmdargs[4])
  outputFolder  <- cmdargs[5]
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnnoTib <- as_tibble(peakAnno@anno)

#combine source data for peaks into one tibble
peakstib <- dplyr::left_join(peaksanno, peaksdat, by = "peakID") %>% 
  dplyr::select(-genomic_feature, -enhType, -H1hesc, -CD56, -CD3, -CD19) %>% 
  dplyr::rename(CD34_dhs = CD34, CD14_dhs = CD14)

peakAnnoTib <- peakAnnoTib %>% 
  dplyr::rename(chr = seqnames) %>% 
  dplyr::mutate(start = start -1)

peakstib <- left_join(peakstib, peakAnnoTib, by = c("chr", "start", "end"))

deseqtib <- deseqtib %>% 
  dplyr::select(symbol, CD34, CD14) %>% 
  dplyr::rename(CD34_rna = CD34, CD14_rna = CD14)

bigtib <- left_join(peakstib, deseqtib) %>% 
  dplyr::mutate(peakFilter = grepl(pattern = "CD[1-3]4", x = accessPattern)) %>% 
  dplyr::filter(peakFilter == TRUE)

#add 1 to avoid infinite values when doing log2 fold changes
bigtib <- bigtib %>% 
  dplyr::mutate(CD34_rna = CD34_rna + 1) %>% 
  dplyr::mutate(CD34_dhs = CD34_dhs + 1) %>% 
  dplyr::mutate(CD14_rna = CD14_rna + 1) %>% 
  dplyr::mutate(CD14_dhs = CD14_dhs + 1)

#annotate max peak value between CD34, CD14 
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(maxValDHS = max(CD14_dhs, CD34_dhs)) %>% 
  dplyr::mutate(maxValRNA = max(CD14_rna, CD34_rna))

bigtib <- bigtib %>% 
  dplyr::filter(maxValRNA > 5)

#adjust annotation formatting
bigtib <- bigtib %>% 
  tidyr::separate(annotation, into = c("anno", NA), sep = "[(]")

bigtib <- bigtib %>% 
  dplyr::mutate(annotation = case_when(anno == "Intron " ~ "gene body",
                                       anno == "Exon " ~ "gene body",
                                       anno == "3' UTR" ~ "gene body",
                                       anno == "5' UTR" ~ "gene body",
                                       anno == "Promoter " ~ "promoter",
                                       anno == "Downstream " ~ "downstream",
                                       TRUE ~ "intergenic"))



plotTib <- bigtib %>%  
  dplyr::select(symbol, annotation, CD34_dhs, CD34_rna, CD14_dhs, CD14_rna)

#remove peaks from genes where the rna value for both samples is zero
plotTib <- plotTib %>% 
  dplyr::mutate(expressionFilter = CD14_rna == 0 & CD34_rna == 0)

plotTib <- plotTib %>% 
  dplyr::filter(expressionFilter  == FALSE)
#replace zero values with 1 for calculating fold changes
plotTib[plotTib == 0] = 1


plotTib <- plotTib %>% 
  dplyr::mutate(log2FC_expression = log2(CD14_rna/CD34_rna)) %>% 
  dplyr::mutate(log2FC_access = log2(CD14_dhs/CD34_dhs))

#find median per anno  per gene
plotTib <- plotTib %>% 
  group_by(symbol, annotation) %>% 
  dplyr::mutate(median_log2FC_access = median(log2FC_access))

plotTib <- plotTib %>% 
  dplyr::select(annotation, median_log2FC_access, log2FC_expression) %>% unique()

p1  <- ggscatter(plotTib, x = "median_log2FC_access", y = "log2FC_expression", facet.by = "annotation", alpha = 0.4, shape = 16, size = 1, add="reg.line", cor.method = "pearson") +
  stat_cor(label.x = -6, label.y = 10)  +
  xlab("log2 FC median Accessibility (norm. fragment counts)") +
  ylab("log2 FC Expression ")

outputLoc <- paste0(outputFolder,"/gonzalez_peakAnno_median.svg")
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

