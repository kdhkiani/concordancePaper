theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  geneTib        <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv')) #karun's input
  outputFolder   <- here("plots", "validation")
} else {
  geneTib        <- read_tsv(cmdargs[1])
  outputFolder   <- cmdargs[2]
}

subsetTib <- geneTib %>% dplyr::select(matches("(_tpm)|(gene_name)")) %>% dplyr::select(- matches("(-and-)|(halfDensity)|(highDensity)")) %>% 
  dplyr::select(- `46-EtOH-nlDensity_tpm`) #sample 46-EtOH-nlDensity is a technical replicate

##For RA

gene <- "HOXA1"

plotTib <- subsetTib %>% dplyr::filter(gene_name == gene) %>% dplyr::select(matches("(gene_name)|(RA)|(EtOH)"))
plotTib <- melt(plotTib) %>% tidyr::separate(col="variable", sep = "-", into = c(NA, "Signal", "Dose")) %>% 
  tidyr::separate(col="Dose", sep="_", into=c("Dose", NA))

plotTib <- plotTib %>% group_by(Dose) %>% dplyr::summarise_at("value", c(mean,sd), na.rm = TRUE) %>% 
  dplyr::rename(mean = fn1, sd = fn2)

plotTib <- plotTib %>% 
  dplyr::mutate(sem = sd/sqrt(3))

p1 <- plotTib %>% 
  dplyr::mutate(Dose = factor(Dose, levels = c("nlDensity", "low", "med", "high"))) %>% 
  ggplot(aes(x=Dose, y=mean, fill=Dose)) +
  geom_bar(stat="identity") +
  geom_linerange(aes(x=Dose, y=mean, ymin=mean-sem, ymax=mean+sem)) +
  scale_fill_manual(values = c("#808782", "#EDABAD", "#DF686C", "#C1292E")) +
  ylab("Mean Expression (TPM)") + xlab("") + theme(legend.position = "none")

outputLoc <- paste0(outputFolder,"/","HOXA1_TPMbarplot.svg")
  
ggsave(filename = outputLoc, plot = p1)



gene <- "SLC5A5"

plotTib <- subsetTib %>% dplyr::filter(gene_name == gene) %>% dplyr::select(matches("(gene_name)|(RA)|(EtOH)"))
plotTib <- melt(plotTib) %>% tidyr::separate(col="variable", sep = "-", into = c(NA, "Signal", "Dose")) %>% 
  tidyr::separate(col="Dose", sep="_", into=c("Dose", NA))

plotTib <- plotTib %>% group_by(Dose) %>% dplyr::summarise_at("value", c(mean,sd), na.rm = TRUE) %>% 
  dplyr::rename(mean = fn1, sd = fn2)

plotTib <- plotTib %>% 
  dplyr::mutate(sem = sd/sqrt(3))

p2 <- plotTib %>% 
  dplyr::mutate(Dose = factor(Dose, levels = c("nlDensity", "low", "med", "high"))) %>% 
  ggplot(aes(x=Dose, y=mean, fill=Dose)) +
  geom_bar(stat="identity") +
  geom_linerange(aes(x=Dose, y=mean, ymin=mean-sem, ymax=mean+sem)) +
  scale_fill_manual(values = c("#808782", "#EDABAD", "#DF686C", "#C1292E")) +
  ylab("Mean Expression (TPM)") + xlab("") + theme(legend.position = "none")

outputLoc <- paste0(outputFolder,"/","SLC5A5_TPMbarplot.svg")

ggsave(filename = outputLoc, plot = p2)

