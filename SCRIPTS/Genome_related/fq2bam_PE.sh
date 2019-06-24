#!/bin/bash

## from Eevi
# Example:
   #####for i in ./raw_data/*_R1_001.fastq.gz; do t=`basename $i _R1_001.fastq.gz`; nice -n 19 ./run_trimgalore_bwa_rmdup.sh $i $t; done
   #####ls ./*_R1_001.fastq.gz| parallel -j 8 t= basename {} _R1_001.fastq.gz\; nice -n 19 fq2bam_PE.sh {} $t
#   ls ./testSet/*_R1_001.fastq.gz | parallel -j 8 nice -n 19 fq2bam_PE.sh {} testout/
[ "$#" -eq "3" ] || { echo "Usage:
ls ./testSet/*_R1_001.fastq.gz | parallel -j 8 nice -n 19 fq2bam_PE.sh {} testout/ 6
# $1 is R1 filename, $2 is output dir, $3 is threads per bwa align
"; exit 1; }



name_no_suffix=`basename $1 _R1_001.fastq.gz`
dir=`dirname $1`
outDir=$2
trimmedDir=$outDir/trimmed/
mkdir -p $trimmedDir
mkdir -p $outDir/withDup
genome="/wrk/data/bwa/hg19/ucsc.hg19.fasta"

set -o nounset
set -o errexit
set -o pipefail

set -x
trim_galore -q 30 --paired --stringency 5 --fastqc --gzip -o $trimmedDir $1 ${dir}/${name_no_suffix}_R2_001.fastq.gz
bwa aln -t $3 -q 20 -f ${outDir}/${name_no_suffix}_val_1.sai $genome ${trimmedDir}${name_no_suffix}_R1_001_val_1.fq.gz
bwa aln -t $3 -q 20 -f ${outDir}/${name_no_suffix}_val_2.sai $genome ${trimmedDir}${name_no_suffix}_R2_001_val_2.fq.gz
bwa sampe $genome ${outDir}/${name_no_suffix}_val_1.sai ${outDir}/${name_no_suffix}_val_2.sai ${trimmedDir}${name_no_suffix}_R1_001_val_1.fq.gz ${trimmedDir}${name_no_suffix}_R2_001_val_2.fq.gz | samtools view -q 20 -b -S - > ${outDir}/${name_no_suffix}.bam
samtools sort ${outDir}/${name_no_suffix}.bam ${outDir}/${name_no_suffix}.sorted
bedtools intersect -abam ${outDir}/${name_no_suffix}.sorted.bam -b ~/DefaultFilesforProc/Genome_related/hg19/hg19_autosomes_human.bed | bedtools intersect -abam stdin -b ~/DefaultFilesforProc/Genome_related/hg19/hg19_blacklist_MergedList.bed -v > ${outDir}/${name_no_suffix}.sorted_autosome_noBlist.bam
samtools rmdup ${outDir}/${name_no_suffix}.sorted_autosome_noBlist.bam ${outDir}/${name_no_suffix}.sorted_rmdup_autosome_noBlist.bam

rm ${trimmedDir}${name_no_suffix}_R1_001_val_1.fq.gz ${trimmedDir}${name_no_suffix}_R2_001_val_2.fq.gz
rm ${outDir}/${name_no_suffix}_val_1.sai ${outDir}/${name_no_suffix}_val_2.sai
rm ${outDir}/${name_no_suffix}.bam
rm ${outDir}/${name_no_suffix}.sorted.bam
mv ${outDir}/${name_no_suffix}.sorted_autosome_noBlist.bam $outDir/withDup/

samtools index ${outDir}/${name_no_suffix}.sorted_rmdup_autosome_noBlist.bam
samtools index $outDir/withDup/${name_no_suffix}.sorted_autosome_noBlist.bam

