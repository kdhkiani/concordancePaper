########
#using source data from paper website
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peaksanno   <- read_tsv(here('Gonzalez2015', 'peaksTable.tsv'))
  peaksdat    <- read_tsv(here('Gonzalez2015', 'DNaseCnts.tsv'))
  deseqtib    <- read_tsv(here('Gonzalez2015', 'RNAseqCnts.tsv')) 
  outputFolder <- here('extractedData', 'gonzalezConcordanceTib.tsv')
} else {
  peaksanno   <- read_tsv(cmdargs[1])
  peaksdat    <- read_tsv(cmdargs[2])
  deseqtib    <- read_tsv(cmdargs[3])
  outputFolder   <- cmdargs[4]
}

peaktib <- left_join(peaksanno, peaksdat)

peaktib <- peaktib %>% 
  dplyr::rename(CD34_dhs = CD34, CD14_dhs= CD14) %>% 
  dplyr::select(- H1hesc, - CD56, - CD3, -CD19)

deseqtib <- deseqtib %>% 
  dplyr::rename(CD34_rna = CD34, CD14_rna= CD14) %>% 
  dplyr::select(- H1hesc, - CD56, - CD3, -CD19)

bigtib <- left_join(peaktib, deseqtib, by = "symbol")

bigtib <- bigtib %>% 
  dplyr::mutate(peakFilter = grepl(pattern = "CD[1-3]4", x = accessPattern)) %>% 
  dplyr::filter(peakFilter == TRUE)

##massage and prepare data from bigtib

#add 1 to avoid infinite values when doing log2 fold changes
bigtib <- bigtib %>% 
  dplyr::mutate(CD34_rna = CD34_rna + 1) %>% 
  dplyr::mutate(CD34_dhs = CD34_dhs + 1) %>% 
  dplyr::mutate(CD14_rna = CD14_rna + 1) %>% 
  dplyr::mutate(CD14_dhs = CD14_dhs + 1)

#annotate max peak value between CD34, CD14 
bigtib <- bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(maxValDHS = max(CD14_dhs, CD34_dhs)) %>% 
  dplyr::mutate(maxValRNA = max(CD14_rna, CD34_rna))

#remove genes where max rna between cell populations is less than 5 tpm
bigtib <- bigtib %>% 
  dplyr::filter(maxValRNA > 5)

#calculate log2FC for accessibility and create a filtering column
bigtib <- bigtib %>% 
  dplyr::mutate(log2FCdhs = log2(CD14_dhs/CD34_dhs)) %>% 
  dplyr::mutate(peakCutoff = abs(log2FCdhs) > 2)

#annotate peak change direction
bigtib <- bigtib %>% 
  dplyr::mutate(peakDir = case_when(CD14_dhs > CD34_dhs ~ "UP",
                                    TRUE ~ "DOWN"))
#determine log2FC for expression
bigtib <- bigtib %>% 
  dplyr::mutate(log2FC = log2(CD14_rna/CD34_rna))

#arrange based on log2FC expression values
bigtib <- bigtib %>% 
  dplyr::arrange(desc(log2FC))

#count number of peak combinations
bigtib <- bigtib %>% 
  group_by(symbol) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & peakCutoff == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & peakCutoff == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & peakCutoff == FALSE))

write_tsv(x = bigtib, file = outputFolder)
