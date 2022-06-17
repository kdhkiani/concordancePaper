joinNearbyPeaksToGenes <- function(deseqtib_path, peakstib_path, window.radius, output.file){
  
  deseqtib     <- read_tsv(deseqtib_path)
  peakstib <- read_tsv(peakstib_path)
  
  peakstib <- peakstib %>%
    mutate(peak_chrom = chrom,
           peak_startLoc = startLocs,
           peak_endLoc   = endLocs) %>%
    dplyr::select(-chrom, -startLocs, -endLocs)
  TSS.window.radius <- as.numeric(window.radius)
  ##########
  
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
}