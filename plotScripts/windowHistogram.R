theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) { #use fread to select specific columns to save memory hassles
  diffPeaks_1kbTib      <- fread(here('extractedData', '1KB_mergedSummit_joinedTablePeaksNearGenes.tsv'), 
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks_2kbTib      <- fread(here('extractedData', '2KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks_5kbTib      <- fread(here('extractedData', '5KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks_10kbTib     <- fread(here('extractedData', '10KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks_25kbTib     <- fread(here('extractedData', '25KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks_50kbTib     <- fread(here('extractedData', '50KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks_100kbTib    <- fread(here('extractedData', '100KB_mergedSummit_joinedTablePeaksNearGenes.tsv'),
                                select = c("gene_name", "peak_chrom", "peak_startLoc", "peak_endLoc", "RA-high_isDeGene", "TGFb-high_isDeGene"))
  diffPeaks             <- read_tsv(here('extractedData', 'differentialAtacPeaksMergedSummitWindows.annotated.tsv'))
  outputLoc             <- here('plots','scatter')
} else {
  diffPeaks_1kbTib      <- read_tsv(cmdargs[1])
  diffPeaks_2kbTib      <- read_tsv(cmdargs[2])
  diffPeaks_5kbTib      <- read_tsv(cmdargs[3])
  diffPeaks_10kbTib     <- read_tsv(cmdargs[4])
  diffPeaks_25kbTib     <- read_tsv(cmdargs[5])
  diffPeaks_50kbTib     <- read_tsv(cmdargs[6])
  diffPeaks_100kbTib    <- read_tsv(cmdargs[7])
  diffPeaks             <- read_tsv(cmdargs[8])
  outputLoc             <- cmdargs[9]
}

#convert these jawns to tibbles
diffPeaks_1kbTib      <- as_tibble(diffPeaks_1kbTib)
diffPeaks_2kbTib      <- as_tibble(diffPeaks_2kbTib)
diffPeaks_5kbTib      <- as_tibble(diffPeaks_5kbTib)
diffPeaks_10kbTib     <- as_tibble(diffPeaks_10kbTib)
diffPeaks_25kbTib     <- as_tibble(diffPeaks_25kbTib)
diffPeaks_50kbTib     <- as_tibble(diffPeaks_50kbTib)
diffPeaks_100kbTib    <- as_tibble(diffPeaks_100kbTib)


#### RA-high first ####

#get diffPeaks 
RA.diffPeaks <- diffPeaks %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":")) %>% 
  dplyr::select(peak, `RA-high-isDiffPeak`) %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE)

### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_1kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                       TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

plotTib <- dplyr::bind_cols(tmp, rep(1, dim(tmp)[1]))
colnames(plotTib)[4] <- "windowSize"

### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_2kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 

  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(2, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)

### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_5kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(5, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)


### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_10kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(10, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)





### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_25kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(25, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)







### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_50kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(50, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)






### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_100kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `RA-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`RA-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`RA-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% RA.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)

#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(100, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)

### plot ###

#alter order and factors for aesthetic ease
plotTib <- plotTib %>% 
  dplyr::mutate(windowSize = as.factor(windowSize)) %>% 
  
  dplyr::mutate(geneAnno = ordered(geneAnno, levels =  c("Non-DEG", "DEG")))

p1 <- ggplot(plotTib, aes(x=windowSize, y=numPeaks, fill=geneAnno)) +
  geom_boxplot(outlier.shape = NA, color = "#F2C945") +
  xlab("Window size (KB)") + ylab("Number of differential peaks winthin window of gene TSS") +
  scale_fill_manual(values = c("#808782", "#C1292E")) + labs(fill = "Gene type") +
  ylim(0,15)

ggsave(filename = paste0(outputLoc,"/RA_sigPeakBoxPlot.svg"), plot = p1)




#### TGFb-high now ####

#get diffPeaks 
TGFb.diffPeaks <- diffPeaks %>% 
  dplyr::mutate(peak = paste(chrom, startLocs, endLocs, sep=":")) %>% 
  dplyr::select(peak, `TGFb-high-isDiffPeak`) %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE)

### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_1kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

plotTib <- dplyr::bind_cols(tmp, rep(1, dim(tmp)[1]))
colnames(plotTib)[4] <- "windowSize"

### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_2kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(2, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)

### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_5kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(5, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)





### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_10kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(10, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)





### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_25kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(25, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)







### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_50kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(50, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)






### Find number of differential peaks per gene for DE and non DE genes ###

tmp <- diffPeaks_100kbTib %>% 
  dplyr::select(gene_name, peak_chrom, peak_startLoc, peak_endLoc, `TGFb-high_isDeGene`) %>% 
  dplyr::mutate(peak = paste(peak_chrom, peak_startLoc, peak_endLoc, sep=":")) %>% 
  dplyr::select(-peak_chrom, -peak_startLoc, -peak_endLoc) %>% 
  
  dplyr::mutate(geneAnno = case_when(`TGFb-high_isDeGene` == 1 ~ "DEG",
                                     TRUE ~ "Non-DEG")) %>% 
  dplyr::select(-`TGFb-high_isDeGene`) %>% 
  
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffPeaks$peak) %>% 
  dplyr::filter(isDiffPeak == TRUE)


#get counts
peakCounts <- tmp %>% group_by(gene_name) %>% 
  dplyr::count(gene_name)

tmp <- left_join(tmp, peakCounts, by="gene_name") %>% 
  dplyr::rename(numPeaks = n) %>% 
  dplyr::select(gene_name, geneAnno, numPeaks) %>% 
  unique()

##add column for window size and add to plotting tibble

tmp.plotTib <- dplyr::bind_cols(tmp, rep(100, dim(tmp)[1]))
colnames(tmp.plotTib)[4] <- "windowSize"

plotTib <- dplyr::bind_rows(plotTib, tmp.plotTib)

### plot ###

#alter order and factors for aesthetic ease
plotTib <- plotTib %>% 
  dplyr::mutate(windowSize = as.factor(windowSize)) %>% 
  
  dplyr::mutate(geneAnno = ordered(geneAnno, levels =  c("Non-DEG", "DEG")))

p2 <- ggplot(plotTib, aes(x=windowSize, y=numPeaks, fill=geneAnno)) +
  geom_boxplot(outlier.shape = NA, color = "#F2C945") +
  xlab("Window size (KB)") + ylab("Number of differential peaks winthin window of gene TSS") +
  scale_fill_manual(values = c("#808782", "#091E5F")) + labs(fill = "Gene type") +
  ylim(0,15)

ggsave(filename = paste0(outputLoc,"/TGFb_sigPeakBoxPlot.svg"), plot = p2)
