##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  annotib     <- read_tsv(here('refs', 'hg19geneAnnotation_Biomart.tsv'))
  deseqtib    <- read_tsv(here('Gonzalez2015', 'RNAseqCnts.tsv')) 
  peaksanno   <- read_tsv(here('Gonzalez2015', 'peaksTable.tsv'))
  peaksdat    <- read_tsv(here('Gonzalez2015', 'DNaseCnts.tsv'))
  output.file  <- here('extractedData', '100KBGonzalez_joinedTablePeaksNearGenes.tsv')
} else {
  deseqtib     <- read_tsv(cmdargs[1])
  peakstib <- read_tsv(cmdargs[2])
  output.file  <- cmdargs[3]
}

#combine source data for peaks into one tibble
peakstib <- dplyr::left_join(peaksanno, peaksdat, by = "peakID")


peakstib <- peakstib %>%
  mutate(peak_chrom = chr,
         peak_startLoc = start,
         peak_endLoc   = end,
         H1hesc_dhs = H1hesc,
         CD34_dhs = CD34,
         CD56_dhs = CD56, 
         CD3_dhs = CD3,
         CD14_dhs = CD14,
         CD19_dhs = CD19) %>%
  dplyr::select(-chr, -start, -end, -H1hesc, -CD34, -CD56, -CD3, -CD14, -CD19)

#filter our annotation table from biomart for GENCODE basic annotations
#annotib <- annotib %>%  dplyr::filter(`GENCODE basic annotation` == "GENCODE basic")

#extract relevant data for annotating by joining expression data with annotation
colnames(deseqtib)[1] <- "gene_name"
colnames(annotib)[15] <- "gene_name"

deseqtib <- left_join(deseqtib, annotib, by="gene_name")

deseqtib <- deseqtib %>%
  mutate(TSS_loc = `Transcription start site (TSS)`,
         chrom = paste0("chr",`Chromosome/scaffold name`),
         strand   = Strand,
         ensg  =  `Gene stable ID`,
         H1hesc_rna = H1hesc,
         CD34_rna = CD34,
         CD56_rna = CD56, 
         CD3_rna = CD3,
         CD14_rna = CD14,
         CD19_rna = CD19) %>%
  dplyr::select(-`Chromosome/scaffold name`, -Strand, -`Gene stable ID`, -`Transcription start site (TSS)`, -H1hesc, -CD34, -CD56, -CD3, -CD14, -CD19)

#drop NA values

deseqtib <- deseqtib %>% drop_na(any_of("TSS_loc"))

#deseqtib <- deseqtib %>% group_by(gene_name) %>% summarise(transcript_start = min(`Transcript start (bp)`), transcript_end = max(`Transcript end (bp)`), chrom = chrom, strand = strand, ensg = ensg,
#                                                           H1hesc_rna = H1hesc_rna, CD34_rna = CD34_rna, CD56_rna = CD56_rna, CD3_rna = CD3_rna, CD14_rna = CD14_rna, CD19_rna = CD19_rna) %>% unique()

TSS.window.radius <- 50000

##########

# grTssWindows <- GRanges(seqnames  = c(deseqtib$chrom),
#                         ranges    = IRanges(start = deseqtib$transcript_start - TSS.window.radius,
#                                             end   = deseqtib$transcript_end + TSS.window.radius),
#                         strand    = deseqtib$strand,
#                         gene_id   = deseqtib$ensg)

grTssWindows <- GRanges(seqnames  = c(deseqtib$chrom),
                        ranges    = IRanges(start = deseqtib$TSS_loc - TSS.window.radius,
                                            end   = deseqtib$TSS_loc + TSS.window.radius),
                        strand    = deseqtib$strand,
                        gene_id   = deseqtib$ensg)

grDiffPeaks  <- GRanges(seqnames = peakstib$peak_chrom,
                        ranges = IRanges(start = peakstib$peak_startLoc,
                                         end   = peakstib$peak_endLoc))

TssWindowPeakOverlaps <- findOverlaps(grTssWindows, grDiffPeaks, ignore.strand=TRUE)
gene_indices <- queryHits(TssWindowPeakOverlaps)
peak_indices <- subjectHits(TssWindowPeakOverlaps)

joinedTibble <- as_tibble(cbind(deseqtib[gene_indices,], peakstib[peak_indices,]))

write_tsv(joinedTibble, output.file, col_names = T)





