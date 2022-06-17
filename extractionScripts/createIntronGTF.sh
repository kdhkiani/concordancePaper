#!/bin/bash

# This is the ource annotation file
# Note that the annotation file should be in Ensembl format (downloaded from Ensembl FTP)
# Specifically, the chromosome names should lack the "chr" prefix
# Also, the first tag in column 9 must be gene_id, and the third tag must be transcript_id
# Also, the source annotation file must have a transcript_source tage in column 9

input="Homo_sapiens.GRCh38.87.chr.gtf"

# 1. Identify transcripts that are supported by both Ensembl and Havana annotations
# 2. Identify and write constitutive exons, i.e. those that appear in all Ensembl/Havana isoforms of a gene
output="./Homo_sapiens.GRCh38.87.chr.consExons.gtf"
cat $input | grep 'transcript_source "ensembl_havana"' | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { split(i,fields,":"); printf("%s\tensembl_havana\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' > $output

# 1. Identify transcripts that are supported by both Ensembl and Havana annotations
# 2. Find all exons
# 3. Write the coordinate of intronic regions, i.e. the regions that separate two adjacent exons of the same gene
output="./Homo_sapiens.GRCh38.87.chr.Introns.gtf"
cat $input | grep 'transcript_source "ensembl_havana"' | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) { split(i,fields,":"); printf("%s\tensembl_havana\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' | bedtools sort -i stdin | awk -v FS='\t' '{ if( last_exon[$9]==1 && (last_exon_end[$9]+1)<($4-1) ) printf("%s\t%s\tintron\t%i\t%i\t%s\t%s\t%s\t%s\n",$1,$2,last_exon_end[$9]+1,$4-1,$6,$7,$8,$9); last_exon[$9]=1; last_exon_end[$9]=$5; }' > $output
