#!/bin/bash

#run with bedtools v.2.15.0 on Penn HPC

dir="/Volumes/TISIPHONE/swami/paper/Analysis_SI2-SI4/extractedData"
inputFile=${dir}/allCondsMergedSummitWindows.bed

bedtools merge -d 50 -i ${inputFile} > ${dir}/allCondsMergedSummitWindows_mergedist50.bed

bedtools merge -d 100 -i ${inputFile} > ${dir}/allCondsMergedSummitWindows_mergedist100.bed

bedtools merge -d 250 -i ${inputFile} > ${dir}/allCondsMergedSummitWindows_mergedist250.bed
