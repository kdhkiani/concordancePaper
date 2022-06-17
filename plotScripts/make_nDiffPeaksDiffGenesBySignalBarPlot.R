cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  deseqTibAnno     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
  peakTibAnno <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows.annotated.tsv'))
  output.file  <- here('plots')
} else {
  deseqtib     <- read_tsv(cmdargs[1])
  peakstib <- read_tsv(cmdargs[2])
  output.file  <- cmdargs[3]
}

ndeg.ralow <- sum(deseqTibAnno$`RA-low_isDeGene`)
ndeg.ramed <- sum(deseqTibAnno$`RA-med_isDeGene`)
ndeg.rahigh <- sum(deseqTibAnno$`RA-high_isDeGene`)
ndeg.tgfblow <- sum(deseqTibAnno$`TGFb-low_isDeGene`)
ndeg.tgfbmed <- sum(deseqTibAnno$`TGFb-med_isDeGene`)
ndeg.tgfbhigh <- sum(deseqTibAnno$`TGFb-high_isDeGene`)

ndeg.vec <- c(ndeg.ralow, ndeg.ramed, ndeg.rahigh, 
              ndeg.tgfblow, ndeg.tgfbmed, ndeg.tgfbhigh)

ord.factor <- c("RA, 50 nM", "RA, 200 nM", "RA, 400 nM", 
                "TGFb, 1.25 ng/mL", "TGFb, 5 ng/mL", "TGFb, 10ng/mL")
barplot.xlabel <- factor(ord.factor, levels = ord.factor)

bptib <- tibble(signal = barplot.xlabel, nDeGenes = ndeg.vec)
p1 <- ggplot(bptib, aes(x = barplot.xlabel, y = nDeGenes)) + 
  geom_bar(stat = "identity", width = .5, color = "black") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4))
p1 <- p1 + xlab("") + ylab("")
p1
ggsave(here('plots', 'validation', 'nDeGenesBySignalPlot.svg'), plot = p1, width = 4, height = 3)

ndiffpeaks.ralow <- sum(peakTibAnno$`RA-low-isDiffPeak`)
ndiffpeaks.ramed <- sum(peakTibAnno$`RA-med-isDiffPeak`)
ndiffpeaks.rahigh <- sum(peakTibAnno$`RA-high-isDiffPeak`)
ndiffpeaks.tgfblow <- sum(peakTibAnno$`TGFb-low-isDiffPeak`)
ndiffpeaks.tgfbmed <- sum(peakTibAnno$`TGFb-med-isDiffPeak`)
ndiffpeaks.tgfbhigh <- sum(peakTibAnno$`TGFb-high-isDiffPeak`)

ndiffpeaks.vec <- c(ndiffpeaks.ralow, ndiffpeaks.ramed, ndiffpeaks.rahigh, 
                    ndiffpeaks.tgfblow, ndiffpeaks.tgfbmed, ndiffpeaks.tgfbhigh)
ord.factor <- c("RA, 50 nM", "RA, 200 nM", "RA, 400 nM", 
                "TGFb, 1.25 ng/mL", "TGFb, 5 ng/mL", "TGFb, 10ng/mL")
barplot.xlabel <- factor(ord.factor, levels = ord.factor)

bptib <- tibble(signal = barplot.xlabel, nDiffPeaks = ndiffpeaks.vec)
p2 <- ggplot(bptib, aes(x = barplot.xlabel, y = nDiffPeaks)) + 
  geom_bar(stat = "identity", width = .5, color = "black") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4))
p2 <- p2 + xlab("") + ylab("")
p2
ggsave(here('plots', 'validation', 'nDiffPeaksBySignalPlot.svg'), plot = p2, width = 4, height = 3)