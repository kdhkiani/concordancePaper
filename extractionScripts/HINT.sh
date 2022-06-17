#!/bin/bash

#run with HINT from the regulatory genomics toolbox

rgt-hint footprinting --atac-seq --paired-end --organism=hg38 --output-location=./HINTatac --output-prefix=EtOH EtOH_merged.bam EtOH.bed

rgt-hint footprinting --atac-seq --paired-end --organism=hg38 --output-location=./HINTatac --output-prefix=RAhigh RAhigh_merged.bam RAhigh.bed

rgt-hint footprinting --atac-seq --paired-end --organism=hg38 --output-location=./HINTatac --output-prefix=TGFb TGFbhigh_merged.bam TGFbhigh.bed

rgt-motifanalysis matching --organism=hg38 --input-files RAhigh.bed EtOH.bed --output-location ./match_RA/

rgt-motifanalysis matching --organism=hg38 --input-files TGFb.bed EtOH.bed --output-location ./match_TGFb

rgt-hint differential --organism=hg38 --bc --nc 30 --mpbs-files=./match_RA/RAhigh_mpbs.bed,./match_RA/EtOH_mpbs.bed --reads-files=../RAhigh_merged.bam,../EtOH_merged.bam --conditions=RAhigh,EtOH --output-location=RAhigh_EtOH

rgt-hint differential --organism=hg38 --bc --nc 30 --mpbs-files=./match_TGFb/TGFb_mpbs.bed,./match_TGFb/EtOH_mpbs.bed --reads-files=../TGFbhigh_merged.bam,../EtOH_merged.bam --conditions=TGFb,EtOH --output-location=TGFb_EtOH
