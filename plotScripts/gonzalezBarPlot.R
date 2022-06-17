theme_set(theme_classic())

########
#using source data from paper website
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksanno   <- read_tsv(here('extractedData', 'Gonzalez2015', 'peaksTable.tsv'))
  peaksdat    <- read_tsv(here('extractedData', 'Gonzalez2015', 'DNaseCnts.tsv'))
  deseqtib    <- read_tsv(here('extractedData', 'Gonzalez2015', 'RNAseqCnts.tsv')) 
  outputLoc <- here('plots', 'validation')
} else {
  peaksanno   <- read_tsv(cmdargs[1])
  peaksdat    <- read_tsv(cmdargs[2])
  deseqtib    <- read_tsv(cmdargs[3])
  outputLoc   <- cmdargs[4]
}

peaktib <- left_join(peaksanno, peaksdat)

peaktib <- peaktib %>% 
  dplyr::rename(CD34_dhs = CD34, CD14_dhs= CD14) %>% 
  dplyr::select(- H1hesc, - CD56, - CD3, -CD19)

deseqtib <- deseqtib %>% 
  dplyr::rename(CD34_rna = CD34, CD14_rna= CD14) %>% 
  dplyr::select(- H1hesc, - CD56, - CD3, -CD19)

bigtib <- left_join(peaktib, deseqtib, by = "symbol")
##massage and prepare data from bigtib

#calculate log2FC for tpm values between CD34,CD14 and create a filtering column
bigtib <- bigtib %>% 
  dplyr::mutate(log2FC = log2(CD14_rna/CD34_rna)) %>% 
  dplyr::mutate(geneCutoff = abs(log2FC) > 2)

#calculate log2FC for accessibility and create a filtering column
bigtib <- bigtib %>% 
  dplyr::mutate(log2FCdhs = log2(CD14_dhs/CD34_dhs)) %>% 
  dplyr::mutate(peakCutoff = abs(log2FCdhs) > 2)

#annotate max peak value between CD34, CD14 
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(maxVal = max(CD14_dhs, CD34_dhs))


#annotate peak change direction
bigtib <- bigtib %>% 
  dplyr::mutate(peakDir = case_when(CD14_dhs > CD34_dhs ~ "Up in CD14",
                                    TRUE ~ "Up in CD34"))
#annotate gene change direction
bigtib <- bigtib %>% 
  dplyr::mutate(geneDir = case_when(CD14_rna > CD34_rna ~ "Up in CD14",
                                    TRUE ~ "Up in CD34"))

#barplot of expression changes
genePlotTib <- bigtib %>% 
  dplyr::filter(geneCutoff == TRUE) %>% 
  dplyr::select(symbol, geneDir) %>% 
  unique() %>% 
  dplyr::count(geneDir)

p1 <- genePlotTib %>% 
  dplyr::mutate(geneDir = factor(geneDir, levels = c("Up in CD34", "Up in CD14"))) %>% 
  ggplot(aes(x=geneDir, y=n, fill=geneDir)) +
  geom_bar(stat="identity") +
  xlab("") + ylab("") + 
  scale_fill_manual(values = c("#3358FF", "#EE8833"))


ggsave(filename = paste0(outputLoc, "/gonzalezExpressionBarPlot.svg"), plot = p1)

#barplot of accessibility changes
peakPlotTib <- bigtib %>% 
  dplyr::filter(peakCutoff == TRUE) %>% 
  dplyr::select(symbol, peakDir) %>% 
  unique() %>% 
  dplyr::count(peakDir)

p2 <- peakPlotTib %>% 
  dplyr::mutate(peakDir = factor(peakDir, levels = c("Up in CD34", "Up in CD14"))) %>% 
  ggplot(aes(x=peakDir, y=n, fill=peakDir)) +
  geom_bar(stat="identity") +
  xlab("") + ylab("") + 
  scale_fill_manual(values = c("#3358FF", "#EE8833"))

ggsave(filename = paste0(outputLoc, "/gonzalezAccessBarPlot.svg"), plot = p2)

