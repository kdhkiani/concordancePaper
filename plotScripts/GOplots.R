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

combined.msigdb <- bind_rows(C5.msigdb, Hallmark.msigdb)

### Retinoic Acid ###

#prepare input
gene_list <- deseqtib %>% 
  dplyr::select(ensg, `RA-high_log2fc`, `RA-high_padj`) %>% 
  dplyr::filter(`RA-high_log2fc` > 2 & `RA-high_padj` < 0.05)

gene_list <- gene_list %>% 
  dplyr::arrange(desc(`RA-high_log2fc`))

#remove NA values and sort in decreasing order
gene_list <- na.omit(gene_list)

genes <- gene_list$ensg

#create gene set enrichment object
enrich <- enricher(gene = genes,
         TERM2GENE = combined.msigdb,
         universe = deseqtib$ensg,
         minGSSize = 3,
         maxGSSize = 1000,
         pvalueCutoff = 0.05, 
         pAdjustMethod = "fdr")

#extract necessary values for plotting manually instead of using out of the box plotting
plotTib <- dplyr::bind_cols(enrich$ID, enrich$p.adjust)
colnames(plotTib) <- c("geneSet", "p.adjust")

plotTib <- plotTib %>% 
  dplyr::mutate(logP = -log10(p.adjust)) %>%  #-log10 transform adjusted p values
  dplyr::select(geneSet, logP) %>% #select relevant columns
  dplyr::slice_max(n = 10, order_by = logP) %>% #select the top ten most signficant
  dplyr::mutate_if(is.character, str_replace_all, pattern="_", replacement=" ") #get rid of '_' in names

p1 <- plotTib %>% 
      dplyr::mutate(geneSet=factor(geneSet, levels=rev(geneSet))) %>% 
      ggplot(aes(x = geneSet, y =logP)) +
      geom_bar(stat = "identity", fill="#C1292E", width = 0.6) +
      coord_flip() + xlab("") + ylab("-log10(p adjusted)")

ggsave(filename = paste0(outputLoc,  "/RAgo.svg"), plot = p1)

### TGF-beta ###

#prepare input
gene_list <- deseqtib %>% 
  dplyr::select(ensg, `TGFb-high_log2fc`, `TGFb-high_padj`) %>% 
  dplyr::filter(`TGFb-high_log2fc` > 2 & `TGFb-high_padj` < 0.05)

gene_list <- gene_list %>% 
  dplyr::arrange(desc(`TGFb-high_log2fc`))

#remove NA values and sort in decreasing order
gene_list <- na.omit(gene_list)

genes <- gene_list$ensg

#create gene set enrichment object
enrich <- enricher(gene = genes,
                   TERM2GENE = combined.msigdb,
                   universe = deseqtib$ensg,
                   minGSSize = 3,
                   maxGSSize = 1000,
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = "fdr")

#extract necessary values for plotting manually instead of using out of the box plotting
plotTib <- dplyr::bind_cols(enrich$ID, enrich$p.adjust)
colnames(plotTib) <- c("geneSet", "p.adjust")

plotTib <- plotTib %>% 
  dplyr::mutate(logP = -log10(p.adjust)) %>%  #-log10 transform adjusted p values
  dplyr::select(geneSet, logP) %>% #select relevant columns
  dplyr::slice_max(n = 10, order_by = logP) %>% #select the top ten most signficant
  dplyr::mutate_if(is.character, str_replace_all, pattern="_", replacement=" ") #get rid of '_' in names

p2 <- plotTib %>% 
  dplyr::mutate(geneSet=factor(geneSet, levels=rev(geneSet))) %>% 
  ggplot(aes(x = geneSet, y =logP)) +
  geom_bar(stat = "identity", fill="#091E5f", width = 0.6) +
  coord_flip() + xlab("") + ylab("-log10(p adjusted)")


ggsave(filename = paste0(outputLoc,  "/TGFbgo.svg"), plot = p2)
