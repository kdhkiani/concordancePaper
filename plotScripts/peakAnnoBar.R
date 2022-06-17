theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  CD14          <- here('extractedData', 'CD14.bed')
  CD34          <- here('extractedData', 'CD34.bed')
  etOHnlDensity <- here('extractedData', 'EtOHnlDensity_atacFragmentCounts_backgroundMinCovg30_mergedist50.bed')
  RA.high       <- here('extractedData', 'RAhigh_atacFragmenCounts_backgroundMinCovg30_mergedist50.bed')
  TGFb.high     <- here('extractedData', 'TGFbhigh_atacFragmentCounts_backgroundMinCovg30_mergedist50.bed')
  outputLoc     <- here('plots', 'validation')
} else {
  CD14          <- cmdargs[1]
  CD34          <- cmdargs[2]
  etOHnlDensity <- cmdargs[3]
  RA.high       <- cmdargs[4]
  TGFb.high     <- cmdargs[5]
  outputLoc <- cmdargs[12]
}

txdb38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

##Gonzalez data 
Gonzalez_peakAnnoList <- lapply(c(CD34, CD14), annotatePeak, TxDb = txdb19,
                                tssRegion = c(-3000,3000), verbose = FALSE)

p1 <- plotAnnoBar(Gonzalez_peakAnnoList)

ggsave(plot = p1, filename = paste0(outputLoc, "/gonzalez_annoBar.svg"))

##Sanford data 
Sanford_peakAnnoList <- lapply(c(etOHnlDensity, RA.high, TGFb.high), annotatePeak, TxDb = txdb38,
                                tssRegion = c(-3000,3000), verbose = FALSE)

p2 <- plotAnnoBar(Sanford_peakAnnoList)

ggsave(plot = p2, filename = paste0(outputLoc, "/sanford_annoBar.svg"))
