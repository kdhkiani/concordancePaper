theme_set(theme_classic())

##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  deseqtib     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
  outputLoc <- here('plots', 'validation')
} else {
  deseqtib     <- read_tsv(cmdargs[1])
  outputLoc  <- cmdargs[2]
}

#load in gene sets from msigdb
C5.msigdb <- msigdbr(species = "Homo sapiens",
                     category = "C5") %>% 
  dplyr::select(gs_name, ensembl_gene)

Hallmark.msigdb <- msigdbr(species = "Homo sapiens",
                           category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)

### Retinoic Acid ###

#prepare input
gene_list <- deseqtib %>% 
  dplyr::select(`RA-high_log2fc`)

gene_list <- gene_list$`RA-high_log2fc`

names(gene_list) <- deseqtib$ensg

#remove NA values and sort in decreasing order
gene_list <- na.omit(gene_list)

gene_list <- sort(gene_list, decreasing = TRUE)

#create gene set enrichment object
gse <- GSEA(geneList = gene_list,
             TERM2GENE = C5.msigdb,
             minGSSize = 3,
             maxGSSize = 1000,
             pvalueCutoff = 0.05,
             verbose = FALSE,
             pAdjustMethod = "fdr")

p1 <- gseaplot2(gse, title = gse$Description[1], geneSetID = 1)

ggsave(filename = paste0(outputLoc,  "/RAgsea.svg"), plot = p1)

###Find p-val and NES to put into svg file in illustrator
gse$NES[1]

gse$p.adjust[1]

### TGF-beta ###

#prepare input
gene_list <- deseqtib %>% 
  dplyr::select(`TGFb-high_log2fc`)

gene_list <- gene_list$`TGFb-high_log2fc`

names(gene_list) <- deseqtib$ensg

#remove NA values and sort in decreasing order
gene_list <- na.omit(gene_list)

gene_list <- sort(gene_list, decreasing = TRUE)

#create gene set enrichment object
gse <- GSEA(geneList = gene_list,
            TERM2GENE = Hallmark.msigdb,
            minGSSize = 3,
            maxGSSize = 1000,
            pvalueCutoff = 0.05,
            verbose = FALSE,
            pAdjustMethod = "fdr")

p2 <- gseaplot2(gse, title = gse$Description[4], geneSetID = 4)

ggsave(filename = paste0(outputLoc,  "/TGFbgsea.svg"), plot = p2)

###Find p-val and NES to put into svg file in illustrator
gse$NES[4]

gse$p.adjust[4]
