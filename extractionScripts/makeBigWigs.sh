#!/bin/bash

#run with samtools v1.1 and deepTools v3.5.1 on Penn HPC


#combine bam files
samtools merge cd14_dnase_combined.final.bam \
./cd14_dnase_r1/cd14_dnase_r1.final.bam \
./cd14_dnase_r2/cd14_dnase_r2.final.bam \
./cd14_dnase_r3/cd14_dnase_r3.final.bam

samtools merge cd34_dnase_combined.final.bam \
./cd34_dnase_r1/cd34_dnase_r1.final.bam \
./cd34_dnase_r2/cd34_dnase_r2.final.bam \
./cd34_dnase_r3/cd34_dnase_r3.final.bam

#index cobmined bam files
samtools index cd14_dnase_combined.final.bam

samtools index cd34_dnase_combined.final.bam

#create bigWigs with counts per million normalization
bamCoverage --normalizeUsing CPM -b cd14_dnase_combined.final.bam \
-o cd14_dnase_combined_CPM.final.bigWig

bamCoverage --normalizeUsing CPM -b cd34_dnase_combined.final.bam \
-o cd34_dnase_combined_CPM.final.bigWig
