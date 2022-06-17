theme_set(theme_classic())

##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib   <- read_tsv(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'))
  diffpeaks <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.tsv'))
  outputFolder <- here('plots', 'summary', '')
} else {
  bigtib    <- read_tsv(cmdargs[1])
  diffpeaks <- read_tsv(cmdargs[2])
  outputFolder  <- cmdargs[3]
}

diffpeaks <- diffpeaks %>% 
  dplyr::mutate(peak =  paste(chrom, startLocs, endLocs, sep=":"))

##count number of peaks for each gene and plot distribution

plotTib <- bigtib %>% dplyr::count(ensg)

p1 <- plotTib %>% 
  ggplot(aes(x=n)) +
  geom_histogram(bins = 30) +
  scale_fill_manual(values = c("808782")) +
  xlim(0,200) + ylab("number of genes") +
  xlab("peaks per gene (using nearest method)")
  
ggsave(plot = p1, filename = paste0(outputFolder,"peaksPerGeneHist.svg"), dpi = 300)

##measure number of peaks differentially up and down in high dose RA and TGFb

tmp <- diffpeaks %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE)

RA.up <- dim(tmp %>% dplyr::filter(`RA-high-avgFoldchange` > 1))[1]
RA.down <-  dim(tmp %>% dplyr::filter(`RA-high-avgFoldchange` < 1))[1]

tmp <- tmp <- diffpeaks %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE)

TGFb.up <- dim(tmp %>% dplyr::filter(`TGFb-high-avgFoldchange` > 1))[1]
TGFb.down <-  dim(tmp %>% dplyr::filter(`TGFb-high-avgFoldchange` < 1))[1]


plotTib <- tibble(cond= c("RA", "RA", "TGFb", "TGFb"),
                  val=c(RA.up, RA.down*-1, TGFb.up, TGFb.down*-1),
                  dir=c("up", "down", "up", "down"))

plotTib <- plotTib %>% 
  dplyr::mutate(id=paste0(cond,"_",dir)) %>% 
  dplyr::mutate(id=factor(id, levels=rev(c("RA_up", "TGFb_up", "RA_down", "TGFb_down"))))
  

p2 <- plotTib %>% 
  ggplot(aes(x=id, y=val, fill=cond)) +
  geom_bar(stat ="identity", width=0.6)+ ylim(-10000, 10000) +
  geom_hline(yintercept=0) + xlab("") + ylab("number of differential peaks") +
  coord_flip() + scale_fill_manual(values = c("#C1292E", "#091E5F"))

ggsave(plot = p2, filename = paste0(outputFolder,"peakChangeDir.svg"), dpi = 300)

##Find the proportion of concordant changes for DEGs

#subset on RA high DEGs

RA.genes <- bigtib %>% 
  dplyr::filter(`RA-high_isDeGene` == TRUE)

RA.genes <- left_join(RA.genes, diffpeaks)

RA.genes <- RA.genes %>% 
  dplyr::select(ensg, `RA-high_log2fc`, `RA-high-isDiffPeak`, `RA-high-avgFoldchange`) %>% 
  replace(is.na(.), FALSE) #this is due to peaks that are not differentially accessible in any condition will not be found in the diffpeaks object and return "NA" during the join
  
RA.genes <- RA.genes %>% 
  dplyr::mutate(gene_incr = `RA-high_log2fc` > 1) %>% 
  dplyr::mutate(peak_incr = `RA-high-avgFoldchange` > 1) %>% 
  dplyr::mutate(same_direction = gene_incr == peak_incr)

RA.genes <- RA.genes %>% 
   group_by(ensg) %>% 
   dplyr::mutate(concordant= case_when(
     (`RA-high-isDiffPeak` == FALSE) ~ TRUE, 
     TRUE ~ same_direction)) %>%
   dplyr::mutate(anypeak = all(concordant == TRUE)) %>% 
   ungroup() %>% 
  dplyr::select(ensg, anypeak) %>% 
  unique()

tmp <- RA.genes %>% 
  dplyr::count(anypeak)

plotTib <- tmp

#subset on TGFb high DEGs

TGFb.genes <- bigtib %>% 
  dplyr::filter(`TGFb-high_isDeGene` == TRUE)

TGFb.genes <- left_join(TGFb.genes, diffpeaks)

TGFb.genes <- TGFb.genes %>% 
  dplyr::select(ensg, `TGFb-high_log2fc`, `TGFb-high-isDiffPeak`, `TGFb-high-avgFoldchange`) %>% 
  replace(is.na(.), FALSE) #this is due to peaks that are not differentially accessible in any condition will not be found in the diffpeaks object and return "NA" during the join

TGFb.genes <- TGFb.genes %>% 
  dplyr::mutate(gene_incr = `TGFb-high_log2fc` > 1) %>% 
  dplyr::mutate(peak_incr = `TGFb-high-avgFoldchange` > 1) %>% 
  dplyr::mutate(same_direction = gene_incr == peak_incr)

TGFb.genes <- TGFb.genes %>% 
  group_by(ensg) %>% 
  dplyr::mutate(concordant= case_when(
    (`TGFb-high-isDiffPeak` == FALSE) ~ TRUE, 
    TRUE ~ same_direction)) %>%
  dplyr::mutate(anypeak = all(concordant == TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(ensg, anypeak) %>% 
  unique()

tmp <- TGFb.genes %>% 
  dplyr::count(anypeak)

plotTib <- plotTib %>% 
  bind_rows(tmp) %>% 
  dplyr::mutate(signal = c(rep("RA",2), rep("TGFb",2)))

#plot

p3 <- plotTib %>% 
  ggplot(aes(x=signal, y=n, fill=anypeak)) +
  geom_bar(position="fill",stat ="identity", width=0.6) + 
  ylab("proportion of DEGs with concordant peak changes")


ggsave(plot = p3, filename = paste0(outputFolder,"degConcordanceProp.svg"), dpi = 300)

