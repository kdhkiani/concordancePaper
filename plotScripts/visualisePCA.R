theme_set(theme_classic())

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  rna.countMatrix        <- readRDS(here('extractedData', 'counts.RNA-seq-matrix.min-count-filtered.rds'))
  atac.countMatrix       <- readRDS(here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.rds'))
  sampleMetadata         <- read_tsv(here('sampleMetadata_SI2-SI4.txt'), col_names=TRUE)
  outputLoc              <- here('plots', 'validation')
} else {
  rna.countMatrix        <- read_tsv(cmdargs[1])
  atac.countMatrix       <- read_tsv(cmdargs[2])
  sampleMetadata         <- read_tsv(cmdargs[3])
  outputLoc              <- cmdargs[4]
}

verbotenSample <- c('08-TGFb-and-RA-low', '09-TGFb-and-RA-med', '11-TGFb-and-RA-high', '13-TGFb-and-RA-low', '18-TGFb-and-RA-med', '21-TGFb-and-RA-high',
                    '28-TGFb-and-RA-med', '29-TGFb-and-RA-low', '35-TGFb-and-RA-high', '46-EtOH-nlDensity', '52-TGFb-and-RA-med',
                    '04-EtOH-highDensity', '07-EtOH-halfDensity', '14-EtOH-highDensity', '20-EtOH-halfDensity', '25-EtOH-halfDensity',
                    '34-EtOH-highDensity')  # '46-EtOH-nlDensity' is a technical replicate of sample 15. remove for differential expression analysis as including it could artificially decrease variance
### RNA-seq ###

##Prepare data
rna.countMatrix <- rna.countMatrix[, !colnames(rna.countMatrix) %in% verbotenSample]

relevantMetadata <- dplyr::filter(sampleMetadata, !sampleID %in% verbotenSample) %>% dplyr::filter(sampleID %in% colnames(rna.countMatrix))

##Plot PCA of RNAseq data

rna.dds <- DESeqDataSetFromMatrix(countData = rna.countMatrix,
                              colData = relevantMetadata,
                              design =  ~ replicate + condition)

rna.vsd <- vst(rna.dds, blind = FALSE)

pcaData <- plotPCA(rna.vsd, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size =7) +
  ggtitle("PCA of RNA-seq data") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#808782", #etOH
                                "#C1292E","#EDABAD","#DF686C", #RA 
                                "#091E5F", "#BCE7fD", "#306FD8")) #TGFbeta

ggsave(filename = paste0(outputLoc,"/RNApca.svg"), plot = p1)


### ATAC-seq ###

##Prepare data

#extract counts as matrix from RangedSummarizedExperiment object 
atac.countMatrix <- as.matrix(assay(atac.countMatrix)) 

#change format of column names to match names from metadata file
colnames(atac.countMatrix) <- as_vector(lapply(strsplit(colnames(atac.countMatrix), '-rep'), function(x) x[1])) 

#filter for samples that are NOT verboten
atac.countMatrix <- atac.countMatrix[, !colnames(atac.countMatrix) %in% verbotenSample]

relevantMetadata <- dplyr::filter(sampleMetadata, !sampleID %in% verbotenSample) %>% dplyr::filter(sampleID %in% colnames(atac.countMatrix))

##Plot PCA of ATACseq data
atac.dds <- DESeqDataSetFromMatrix(countData = atac.countMatrix,
                         colData = relevantMetadata,
                         design = ~ replicate + condition)

atac.vsd <- vst(atac.dds, blind = FALSE)

pcaData <- plotPCA(atac.vsd, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p2 <- ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size =7) +
  ggtitle("PCA of ATAC-seq data") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = "NA") + #will use same legend as above in figure
  scale_color_manual(values = c("#808782", #etOH
                                 "#C1292E","#EDABAD","#DF686C", #RA 
                                 "#091E5F", "#BCE7fD", "#306FD8")) #TGFbeta

ggsave(filename = paste0(outputLoc,"/ATACpca.svg"), plot = p2)
