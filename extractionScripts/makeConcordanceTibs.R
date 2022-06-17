##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  bigtib   <- read_tsv(here('extractedData', 'tmp_20211122' , 'genePeak_1to1mapping_mergedist50.tsv'))
  diffpeaks <- read_tsv(here('extractedData', 'tmp_20211122' , 'differentialAtacPeaksMergedSummitWindows_minCovg30_minFc1.5_mergedist50.tsv'))
  outputLoc <- here('extractedData','concordanceTib')
} else {
  bigtib    <- read_tsv(cmdargs[1])
  diffpeaks <- read_tsv(cmdargs[2])
  outputLoc  <- cmdargs[3]
}

##massage and prepare data from bigtib

###################
## RA high dose ##

#pull out the relevant values 
RA.bigtib <- bigtib %>% 
  dplyr::select(gene_name, ensg, peak, `RA-high_isDeGene`,
                `EtOH-nlDensity_avgTPM`, `RA-high_avgTPM`, `EtOH-nlDensity-avgNormFragmentCounts`, `RA-high-avgNormFragmentCounts`,
                `RA-med-avgNormFragmentCounts`, `RA-low-avgNormFragmentCounts`)

#calculate log2FC for tpm values between RA/etOH
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`RA-high_avgTPM` = `RA-high_avgTPM` + 1) %>% 
  dplyr::mutate(log2FC = log2(`RA-high_avgTPM`/`EtOH-nlDensity_avgTPM`))


#annotate max peak value between RA/etOH 
RA.bigtib <- RA.bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(RA.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`RA-high-avgNormFragmentCounts`))

#annotate peak change direction
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(peakDir = case_when(`RA-high-avgNormFragmentCounts` > `EtOH-nlDensity-avgNormFragmentCounts` ~ "UP",
                                    TRUE ~ "DOWN"))

##get the differential peak identifiers for RA-high
RA.diffpeaks <- diffpeaks %>% 
  dplyr::mutate(peak = paste(chrom,startLocs,endLocs,sep=":")) %>% 
  dplyr::select(peak, `RA-high-isDiffPeak`) %>% 
  dplyr::filter(`RA-high-isDiffPeak` == TRUE)

##annotate if peaks in bigtib are diffpeaks
RA.bigtib <- RA.bigtib %>% 
  dplyr::mutate(isDiffPeak = peak %in% RA.diffpeaks$peak)

#filter out max peak val across lows < 20
RA.bigtib <- RA.bigtib %>%
  dplyr::filter(RA.maxVal > 30)

#count number of peak combinations
RA.bigtib <- RA.bigtib %>% 
  group_by(gene_name) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & isDiffPeak == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & isDiffPeak == FALSE))

write_tsv(x = RA.bigtib, file = paste0(outputLoc, "_RAhigh.tsv"))

##########################
## TGFb high dose ##

#pull out the relevant values 
TGFb.bigtib <- bigtib %>% 
  dplyr::select(gene_name, ensg, peak, `TGFb-high_isDeGene`,
                `EtOH-nlDensity_avgTPM`, `TGFb-high_avgTPM`, `EtOH-nlDensity-avgNormFragmentCounts`, `TGFb-high-avgNormFragmentCounts`,
                `TGFb-med-avgNormFragmentCounts`, `TGFb-low-avgNormFragmentCounts`)

#calculate log2FC for tpm values between TGFb/etOH
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(`EtOH-nlDensity_avgTPM` = `EtOH-nlDensity_avgTPM` + 1) %>% 
  dplyr::mutate(`TGFb-high_avgTPM` = `TGFb-high_avgTPM` + 1) %>% 
  dplyr::mutate(log2FC = log2(`TGFb-high_avgTPM`/`EtOH-nlDensity_avgTPM`))


#annotate max peak value between TGFb/etOH 
TGFb.bigtib <- TGFb.bigtib %>% 
  rowwise() %>% 
  dplyr::mutate(TGFb.maxVal = max(`EtOH-nlDensity-avgNormFragmentCounts`,`TGFb-high-avgNormFragmentCounts`))

#annotate peak change direction
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(peakDir = case_when(`TGFb-high-avgNormFragmentCounts` > `EtOH-nlDensity-avgNormFragmentCounts` ~ "UP",
                                    TRUE ~ "DOWN"))

##get the differential peak identifiers for TGFb-high
TGFb.diffpeaks <- diffpeaks %>% 
  dplyr::mutate(peak = paste(chrom,startLocs,endLocs,sep=":")) %>% 
  dplyr::select(peak, `TGFb-high-isDiffPeak`) %>% 
  dplyr::filter(`TGFb-high-isDiffPeak` == TRUE)

##annotate if peaks in bigtib are diffpeaks
TGFb.bigtib <- TGFb.bigtib %>% 
  dplyr::mutate(isDiffPeak = peak %in% TGFb.diffpeaks$peak)

#filter out max peak val across lows < 20
TGFb.bigtib <- TGFb.bigtib %>%
  dplyr::filter(TGFb.maxVal > 30)

#count number of peak combinations
TGFb.bigtib <- TGFb.bigtib %>% 
  group_by(gene_name) %>% 
  dplyr::mutate(Up_Sig = sum(peakDir == "UP" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Up_NotSig = sum(peakDir == "UP" & isDiffPeak == FALSE)) %>% 
  dplyr::mutate(Down_Sig = sum(peakDir == "DOWN" & isDiffPeak == TRUE)) %>% 
  dplyr::mutate(Down_NotSig = sum(peakDir == "DOWN" & isDiffPeak == FALSE))

write_tsv(x = TGFb.bigtib, file = paste0(outputLoc, "_TGFbhigh.tsv"))
