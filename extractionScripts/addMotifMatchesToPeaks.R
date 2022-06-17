#adapted from Sanford et al., 2020 script "addIntegrationMetricsAndMotifMatchesToDiffPeaks.R"

register(MulticoreParam(4, progressbar = TRUE))

data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper

## user-defined parameters
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  diffPeakTib    <- read_tsv(here('extractedData','atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'))
  dpOutputTibLoc <- here('extractedData','atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.annotated.tsv')
  selected.PWM.objects <- read_rds(here('extractedData', 'most_variable_motifs_Robject.rds'))
} else {
  diffPeakTib    <- read_tsv(cmdargs[1])
  selected.PWM.objects <- read_rds(cmdargs[2])
  dpOutputTibLoc <- cmdargs[3]
}


# motif stuff...

## body of script
dpOutputTib <- diffPeakTib


peak_gr <- GRanges(seqnames = dpOutputTib$chrom, 
                       ranges = IRanges(start = dpOutputTib$startLocs,
                                        end   = dpOutputTib$endLocs))

# note: details on match scores are here https://github.com/jhkorhonen/MOODS/wiki/Brief-theoretical-introduction
# "Intuitively, this score compares the probability that the model specified by the original PWM generated the sequence 
# versus the probability that the background model generated the sequence"
# note #2: the score sums the log probability ratio of PWM value over background value for each position in the motif. because of this,
# i think it's important to remember that longer motifs can obtain higher scores than shorter motifs.
selectedMotifsMatches <- matchMotifs(selected.PWM.objects, 
                                     peak_gr, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38, 
                                     bg = "genome",
                                     out = "scores") 
motif.names <- sapply(strsplit(colnames(selectedMotifsMatches), "_"), function (x) x[3])
counter <- 1
for (motif.name in motif.names) {
  dpOutputTib[[paste0(motif.name, '_motifMatchScore')]] <- motifScores(selectedMotifsMatches)[, counter]
  dpOutputTib[[paste0(motif.name, '_numMotifMatches')]] <- motifCounts(selectedMotifsMatches)[, counter]
  counter <- counter + 1
}

# now categorize motifs into sets and add motifsets
# these were the top 72 motifs selected by this script, which i then manually categorized into groups

# factors that go down in both RA and TGFb (i.e. EtOH factors)
klf.factors        <- c("KLF1", "KLF2", "KLF3", "KLF4", "KLF5", "KLF7")
ap2.factors        <- c("TFAP2B", "TFAP2C", "TFAP2E", "TFAP2A")
other.etoh.factors <- c("GRHL1", "SP3")
all.etoh.factors   <- c(klf.factors, ap2.factors, other.etoh.factors)

# NFKB factors, which go up a little in both RA and TGFb treatment
nfkb.factors <- c("REL", "RELA", "NFKB1")

# factors that go up in RA treatment
fox.factors      <- c("FOXJ2", "FOXA2", "FOXA1", "FOXA3", "FOXC2", "FOXD2", "FOXD3")
hox.factors      <- c("HOXC13", "HOXD13", "HOXB13", "HOXC10", "HOXA13")
rar.factors      <- c("RARA")
irf.factors.ra   <- c("IRF2", "IRF3", "IRF8")
elf.factors      <- c("ELF1", "ELF2", "ELF3", "ELF4", "ELF5")
other.ra.factors <- c("BCL11A", "BCL11B", "CDX1", "CDX2", "EHF", "GABPA", "NR2C2", "ETV7", "SPI1", "SPIB", "SPIC", "ELK4", "ETS2", "STAT2")
all.ra.factors   <- c(fox.factors, hox.factors, rar.factors, irf.factors.ra, elf.factors, other.ra.factors)

# factors that go up in TGFb treatment
ap1.factors        <- c("FOS", "FOSL1", "FOSL2", "FOSB", "JUN", "JUNB", "JUND", "JDP2", "BATF")
smad.factors       <- c("SMAD3", "SMAD4", "SMAD9")
irf.factors.tgfb   <- c("IRF4")
other.tgfb.factors <- c("BACH1", "BACH2", "SMARCC1", "NFE2", "NFE2L2", "MAFF", "MAFK")
all.tgfb.factors   <- c(ap1.factors, smad.factors, irf.factors.tgfb, other.tgfb.factors)


ra.dominant.factors   <- c("RARA", 'FOXA1', 'FOXA2', 'FOXA3', 'FOXC2', 'FOXD3', 'SPI', 'SPIB', "SPIC", "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5")
tgfb.dominant.factors <- c("SMAD3", "SMAD4", "SMAD9", "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", "BACH1", "BACH2", "BATF",
                           "SMARCC1", "NFE2", "NFE2L2", "MAFF", "MAFK")


factor.group.list  <- list(klf.factors, ap2.factors, other.etoh.factors, all.etoh.factors, nfkb.factors, fox.factors, hox.factors, rar.factors, irf.factors.ra, elf.factors, other.ra.factors, all.ra.factors, 
                           ap1.factors, smad.factors, irf.factors.tgfb, other.tgfb.factors, all.tgfb.factors,
                           ra.dominant.factors, tgfb.dominant.factors)
factor.group.names <-   c('group-KLF', 'group-AP2', 'group-otherEtOH',  'group-allEtOH',  'group-NFKB', 'group-FOX', 'group-HOX', 'group-RAR', 'group-raIRF',  'group-ELF', 'group-otherRA',  'group-allRA',  
                          'group-AP1', 'group-SMAD', 'group-tgfbIRF',  'group-otherTGFB',  'group-allTGFB',
                          'group-RAdominant', 'group-TGFbdominant')
stopifnot(length(factor.group.list) == length(factor.group.names))

for (ii in 1:length(factor.group.names)) {
  factor.group.name <- factor.group.names[ii]
  factor.group <- factor.group.list[[ii]]
  motif_selection_regex <- paste0("(", paste(factor.group, collapse="|"), ")", "_motifMatchScore")
  motifGroupMatchTib <- dplyr::select(dpOutputTib, matches(motif_selection_regex))
  motifGroupMaxScoresAsVector <- sapply(1:nrow(motifGroupMatchTib), function(x) max(as_vector(motifGroupMatchTib[x, ])))
  dpOutputTib[[paste0(factor.group.name, '_maxMotifMatchScore')]] <- motifGroupMaxScoresAsVector
}

###

write_tsv(dpOutputTib, dpOutputTibLoc)
