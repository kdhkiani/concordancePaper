#!/bin/bash

for i in HiCStein-MCF7*; do
    sed -i '' '1d' $i;
done

cat HiCStein-MCF7* > HiCStein-MCF7_allChr_hg19.insulation.boundaries

intersect -wa -wb -a RAhigh_atacFragmenCounts_backgroundMinCovg30_mergedist50.bed -b HiCStein-MCF7_allChr_hg38.bed > TADkey_RAhigh.tsv

intersect -wa -wb -a TGFbhigh_atacFragmentCounts_backgroundMinCovg30_mergedist50.bed -b HiCStein-MCF7_allChr_hg38.bed > TADkey_TGFbhigh.tsv

bedtools intersect -wa -wb -a ../refs/tss.bed -b HiCStein-MCF7_allChr_hg38.bed > tss_TADkey.tsv
