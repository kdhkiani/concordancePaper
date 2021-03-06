#!/bin/bash
# this if statement is needed for bsub
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
source /project/arjunrajlab/resources/atac_pipeline_virtual_env/bin/activate  # this loads a virtual environment that contains macs2 2.1.1.20160309
module load bowtie2/2.3.4.1
module load java/openjdk-1.8.0
module load samtools-1.1
module load R/3.1.1
module load python/2.7.9
module load picard/1.96
module load FastQC-0.11.2
module load homer-v4.6
PATH=\
/project/arjunrajlab/resources/bedtools2/bin:\
/project/arjunrajlab/resources/bedGraphToBigWig:\
$PATH \
$@
module unload bowtie2/2.3.4.1
module unload java/openjdk-1.8.0
module unload samtools-1.1
module unload R/3.1.1
module unload python/2.7.9
module unload picard/1.96
module unload FastQC-0.11.2
module unload homer-v4.6
deactivate
