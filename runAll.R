## Run all extraction scripts
setwd("/Volumes/ERINYES/atac_paper/Analysis_SI2-SI4/")
library(here)

#navigate to the repo -- this will be need to be changed for another computer

#source some helper functions
source(here('extractionScripts', 'util.R'))

#Load libraries in
loadLibraries()

#make consensus peak files

source(here('extractionScripts', 'makeConsensusPeakFilesForEachCondition.R'))

#merge consensus peak fliles

source(here('extractionScripts', 'mergeConsensusPeaksOfSeveralConditions.R'))

#merge bed files using varying distances by executing shell script -- note this will require bedtools

system(here('extractionScripts/combineBamFiles.sh'))

#merge bam files for use with HINT

system(here('extractionScripts/mergeBedFiles.sh'))

#make atac frag count file for consensus peak file of merge distance = 50

source(here('extractionScripts', 'createAtacFragmentCountMatrix.R'))

#Run empirical Fdr based differential peak selection with minimumn coverage of 30 norm. fragment counts
# and minimum fold change of 1.5

source(here('extractionScripts', 'runEmpiricalFdrBasedDifferentialPeakSelection.R'))

#re-run code with different parameter set: min coverage = 10 norm. fragment counts,

cmdargs <- c("10",
             "1.5",
             here('extractedData','atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.rds'),
             here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg10_minFc1.5_mergedist50.tsv'),
             here('plots', 'differentialAtacPeaks_FDR_plot_minCovg10'))

source(here('extractionScripts', 'runEmpiricalFdrBasedDifferentialPeakSelection.R'))

rm(cmdargs)

#make annotation bed files for merged peaks

source(here('extractionScripts', 'makeMergedPeakAnnoBedFiles.R'))

#make annotation bed files for differenential peaks

source(here('extractionScripts', 'makeDiffPeakAnnoBedFiles.R'))

#make annotation bed files for data from Gonzalez et al., 2015
source(here('extractionScripts', 'makeGonzalezAnnoBedFiles.R'))

#make the list of the most variable motifs in the differential peak set

source(here('extractionScripts', 'makeMostVariableMotifSet.R'))

#annotate the peak set with motif matches

source(here('extractionScripts', 'addMotifMatchesToPeaks.R'))

##create one-to-one mapping tibble for 'nearest approach'

#annotate peaks

source(here('extractionScripts', 'peakAnno.R'))


#combine files into one-to-one mapping

source(here('extractionScripts', 'joinGenesToPeaks_1to1table.R'))

#create joinedTibbles for window approach

source(here('extractionScripts', 'joinNearbyPeaksToGenes.R'))


cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "500",
             here('extractedData', '1KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "1500",
             here('extractedData', '3KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "2500",
             here('extractedData', '5KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "5000",
             here('extractedData', '10KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "12500",
             here('extractedData', '25KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.annotated.tsv'),
             '20000',
             here('extractedData', '40KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "25000",
             here('extractedData', '50KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

cmdargs <- c(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'),
             here('extractedData', 'atacFragmentCountsAllCondsMergedSummitWindows_mergedist50.tsv'),
             "50000",
             here('extractedData', '100KB_mergedist50_joinedTablePeaksNearGenes.tsv'))

joinNearbyPeaksToGenes(cmdargs[1], cmdargs[2], cmdargs[3], cmdargs[4])

rm(cmdargs)

source(here('extractionScripts', 'GonzalezjoinNearbyPeaksToGenes.R'))

##make tibble for only distal peaks

source(here('extractionScripts', 'makeDistalPeakTib.R'))

##create tibble for concordance plots

source(here('extractionScripts', 'makeConcordanceTibs.R'))

source(here('extractionScripts', 'gonzalezMakeConcordanceTib.R'))

##Perform HINT analysis

#make bed file for only peaks that are not differentially accessible

source(here('extractionScripts', 'hint_makeNonDiffPeakBed.R'))

system(here('extractionScripts/mergeBedFiles.sh'))

##Process ramirez HL60 to monocyte-derived macrophage data

#get fragment counts

source(here('extractionScripts', 'ramirez_makeConsensusPeakFilesForEachCondition.R'))

source(here('extractionScripts', 'ramirez_MergeConsensusPeaksOfSeveralConditions.R'))

source(here('extractionScripts', 'ramirez_createAtacFragmentCountMatrix.R'))

source(here('extractionScripts', 'ramirez_NormalizeFragmentCounts.R'))

#annotate peaks

source(here('extractionScripts', 'ramirez_peakAnno.R'))

#combine with expression data using 'nearest' method

source(here('extractionScripts', 'ramirez_joinGenesToPeaks_1to1table.R'))

#create bed file of peaks to use to annotate peak location distribution with peakAnnoBarPlot
source(here('extractionScripts', 'ramirez_joinGenesToPeaks_1to1table.RmakeRamirezMergedPeakAnnoBedFiles.R'))

##Process IL1alpha

#get fragment counts

source(here('extractionScripts', 'jw_makeConsensusPeakFilesForEachCondition.R'))

source(here('extractionScripts', 'jw_MergeConsensusPeaksOfSeveralConditions.R'))

source(here('extractionScripts', 'jw_createAtacFragmentCountMatrix.R'))

source(here('extractionScripts', 'jw_NormalizeFragmentCounts.R'))

#annotate peaks

source(here('extractionScripts', 'jw_peakAnno.R'))

#combine with expression data using 'nearest' method

source(here('extractionScripts', 'jw_joinGenesToPeaks_1to1table.R'))

#add motif match data

source(here('extractionScripts', 'jw_addMotifMatchesToPeaks.R'))

#TAD analysis

system(here('extractionScripts/TADprocessing.sh'))

#Perform intorn/exon analysis 
system(here('extractionScripts/createIntronGTF.sh'))

system(here('extractionScripts/intronHTSeq.sh'))
#################################################

###plots###

##get stats for barplots for differentially expessed and accessible genes

source(here('plotScripts', 'make_nDiffPeaksDiffGenesBySignalBarPlot.R'))

##perform GO analysis on differentially upregulated

source(here('plotScripts', 'GOplots.R'))

##perform GSEA analysis

source(here('plotScripts', 'gseaPlots.R'))

##perform motif analysis

source(here('plotScripts', 'makeMotifAnalysisPlots.R'))

##get numbers for Venn diagrams -- Venn circles made in illustrator

source(here('plotScripts', 'getVennNumbers.R'))

##make expression barplots

source(here('plotScripts', 'makeTPMbarplots.R'))

##visualize PCA of RNAseq and ATAC-seq counts

source(here('plotScripts', 'visualisePCA.R'))

##make accessibility plots

source(here('plotScripts', 'peakAccess_ParCoord.R'))

#make summary plots
source(here('plotScripts', 'makeSummaryPlots.R'))

##make density plot of complexity distributions, peak widths by complexity, and gene expression by complexity

source(here('plotScripts', 'lociComplexityDist.R'))

##make proportion plot for MCF-7 data

source(here('plotScripts', 'sanfordMakeProportionPlot.R'))

##make proportion plot for MCF-7 data with 10 norm fragment cutoff

cmdargs <- c(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'),
             here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg10_minFc1.5_mergedist50.tsv'),
             here('plots'))

source(here('plotScripts', 'sanfordMakeProportionPlot.R'))

rm(cmdargs)

##make proportion plot for hematopoietic differentiation data

source(here('plotScripts', 'gonzalezMakeProportionPlot.R'))

##make proportion plot for HL60 differentiation data

source(here('plotScripts', 'ramirezMakeProportionPlot.R'))
##make concordance plots

source(here('plotScripts', 'concordancePlot.R'))

##make scatterplots

#use window method

source(here('plotScripts', 'scatterPlots.R'))

cmdargs <- c(here('extractedData', '100KB_mergedist50_joinedTablePeaksNearGenes.tsv'),
             here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg10_minFc1.5_mergedist50.tsv'),
             here('plots/100kb_'))

sanfordScatterPlots(cmdargs[1], cmdargs[2], cmdargs[3])

rm(cmdargs)

cmdargs <- c(here('extractedData', '40KB_mergedist50_joinedTablePeaksNearGenes.tsv'),
             here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg10_minFc1.5_mergedist50.tsv'),
             here('plots/distal40kb_'))

sanfordScatterPlots(cmdargs[1], cmdargs[2], cmdargs[3])

rm(cmdargs)

cmdargs <- c(here('extractedData', '3KB_mergedist50_joinedTablePeaksNearGenes.tsv'),
             here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg10_minFc1.5_mergedist50.tsv'),
             here('plots/3kb_'))

sanfordScatterPlots(cmdargs[1], cmdargs[2], cmdargs[3])

rm(cmdargs)

#use nearest method

cmdargs <- c(here('extractedData', 'genePeak_1to1mapping_mergedist50.tsv'),
             here('extractedData', 'differentialAtacPeaksMergedSummitWindows_minCovg10_minFc1.5_mergedist50.tsv'),
             here('plots'))

source(here('plotScripts', 'scatterPlots.R'))

rm(cmdargs)

#scatters by peak annotation

source(here('plotScripts', 'scatterByPeakAnno.R'))

#use window method for gonzalez

cmdargs <- c(here('extractedData', 'Gonzalez2015', 'RNAseqCnts.tsv'),
            here('extractedData', 'Gonzalez2015', 'peaksTable.tsv'),
            here('extractedData', 'Gonzalez2015', 'DNaseCnts.tsv'),
            TRUE,
            here('plots', '100kb_'))

source(here('plotScripts', 'gonzalezScatterPlots.R'))

rm(cmdargs)

#use nearest method

source(here('plotScripts', 'gonzalezScatterPlots.R'))

#scatter by peak annotations

source(here('plotScripts', 'gonzalezScatterByPeakAnno.R'))

#scatter using MCF-7 TAD data

source(here('plotscripts', 'TADscatter.R'))

##make window histogram of number of DEGs and non-DEGs per window size

source(here('plotScripts', 'windowHistogram.R'))

##make bar plots of number of differentially expressed genes and accessible peaks from gonzalez

source(here('plotScripts', 'gonzalezBarPlot.R'))

##make bar plots with proportion of peak annotations
source(here('plotScripts', 'peakAnnoBar.R'))

##make scatterplots by motif

source(here('pootScripts', 'motifScatter.R'))

##make scatterplots by intron/exon

source(here('pootScripts', 'intron_scatter.R'))


##make proportion plot and motif scatter for IL1alpha data

source(here('pootScripts', 'jwMakeProportionPlot.R'))

source(here('pootScripts', 'jw_motifScatter.R'))