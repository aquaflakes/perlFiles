#!/bin/bash

## from Eevi
[ "$#" -eq "3" ] || { echo "Usage:
ls ./testSet/*.fastq.gz | parallel -j 8 nice -n 19 fq2bam_SR.sh {} testout/ 6
# $1 is R1 filename, $2 is output dir, $3 is threads per bwa align
"; exit 1; }



name_no_suffix=`basename $1 _001.fastq.gz`
dir=`dirname $1`
outDir=$2
trimmedDir=$outDir/trimmed/
mkdir -p $trimmedDir
genome="/wrk/data/bwa/hg19/ucsc.hg19.fasta"

set -o nounset
set -o errexit
set -o pipefail

set -x
trim_galore -q 30 --stringency 5 --fastqc --gzip -o $trimmedDir $1 
bwa aln -t $3 -q 20 -f ${outDir}/${name_no_suffix}_trimmed.sai $genome ${trimmedDir}${name_no_suffix}_001_trimmed.fq.gz

bwa samse $genome ${outDir}/${name_no_suffix}_trimmed.sai ${trimmedDir}${name_no_suffix}_001_trimmed.fq.gz | samtools view -q 20 -b -S - > ${outDir}/${name_no_suffix}.bam
samtools sort ${outDir}/${name_no_suffix}.bam ${outDir}/${name_no_suffix}.sorted
samtools rmdup ${outDir}/${name_no_suffix}.sorted.bam ${outDir}/${name_no_suffix}.sorted.rmdup.bam

rm ${trimmedDir}${name_no_suffix}_001_trimmed.fq.gz
rm ${outDir}/${name_no_suffix}_trimmed.sai
rm ${outDir}/${name_no_suffix}.bam
rm ${outDir}/${name_no_suffix}.sorted.bam

bedtools intersect -abam ${outDir}/${name_no_suffix}.sorted.rmdup.bam -b ~/DefaultFilesforProc/Genome_related/hg19/hg19_autosomes_human.bed | bedtools intersect -abam stdin -b ~/DefaultFilesforProc/Genome_related/hg19/hg19_blacklist_MergedList.bed -v > ${outDir}/${name_no_suffix}.sorted_rmdup_autosome_noBlist.bam

rm ${outDir}/${name_no_suffix}.sorted.rmdup.bam
samtools index ${outDir}/${name_no_suffix}.sorted_rmdup_autosome_noBlist.bam