#!/bin/bash

#run with bedtools v.2.15.0 on Penn HPC

#comnbined leslie files into consensus peak set
cat CD14_hg38.bed CD34_hg38.bed > leslie.bed

#sort and merge file using a distance parameter of 50 base pairs
bedtools sort -i leslie.bed > leslie.sorted.bed

bedtools merge -i leslie.sorted.bed -d 50 > leslie_mergedist50.bed

#find the peaks where there is an intersection between sandford and leslie
bedtools intersect -a allCondsMergedSummitWindows_mergedist50.bed -b leslie_mergedist50.bed > intersect.bed

#merge peaks together within 50 bps of each other
bedtools merge -i intersect.bed -d 50 > intersect_mergedist50.bed

#compare the relative number of peaks
wc -l intersect_mergedist50.bed

wc -l leslie_mergedist50.bed
