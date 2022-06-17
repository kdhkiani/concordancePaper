#!/bin/bash

#run with samtootools v1.11 on Penn HPC

samtools merge EtOH_merged.bam \
  05-EtOH-nlDensity-rep1.final.bam \
  15-EtOH-nlDensity-rep2.final.bam \
  27-EtOH-nlDensity-rep3.final.bam

samtools index EtOH_merged.bam

samtools merge RAhigh_merged.bam \
  06-RA-high-rep1.final.bam \
  23-RA-high-rep2.final.bam \
  32-RA-high-rep3.final.bam

samtools index RAhigh_merged.bam

samtools merge TGFbhigh_merged.bam \
  10-TGFb-high-rep1.final.bam \
  22-TGFb-high-rep2.final.bam \
  31-TGFb-high-rep3.final.bam

samtools index TGFbhigh_merged.bam
