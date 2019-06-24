
#!/bin/bash

outDir=.
cores=20

saveNAME=${outDir}/cmd.sh
mkdir -p $outDir
cp ~/Nut_zhuData/lib/SCRIPTS/cmd1.sh $saveNAME # back up cmd


Runfolder=170829_K00110_0138_BHJG3CBBXX
thisProject=FJ4.6CAPFL

mkdir -p ~/bckup_Mybook/fastq/$thisProject/Original_reads/
ln -s /fastq/data/$Runfolder/${Runfolder}_Fangjie ~/bckup_Mybook/fastq/$thisProject/Original_reads/
mkdir -p ~/seqFiles2/$thisProject/
mkdir -p ~/seqFiles2/$thisProject/allreads
mkdir -p ~/seqFiles2/$thisProject/analysis
ln -s ~/bckup_Mybook/fastq/$thisProject/Original_reads ~/seqFiles2/$thisProject/

# pear
cd ~/bckup_Mybook/fastq/$thisProject/Original_reads/
for subdir in $(ls -d ${Runfolder}_Fangjie/* | grep -vE "Reports|Stats" | xargs -n1 basename);\
do mkdir -p peared_$subdir;\
   pearAll.pl -f "${Runfolder}_Fangjie/"$subdir"/Trulig147*.gz" -o peared147_$subdir -keep 101 -t 40 && cp peared147_$subdir/* peared_$subdir/ -r && rm peared147_$subdir/ -r;\
   pearAll.pl -f "${Runfolder}_Fangjie/"$subdir"/Trulig200*.gz" -o peared200_$subdir -keep 154 -t 40 && cp peared200_$subdir/* peared_$subdir/ -r && rm peared200_$subdir/ -r;\
   rm -r peared_$subdir/tmp

done;
	

# trim adaptor
cd /wrk/data/fangjie/seqFiles2/$thisProject/allreads
for subdir in $(ls -d ~/bckup_Mybook/fastq/$thisProject/Original_reads/${Runfolder}_Fangjie/* | grep -vE "Reports|Stats" | xargs -n1 basename);\
do FJ4.4_adaptorTrim.pl -f '../Original_reads/peared_'$subdir'/*.gz' -o Adpt_Trimmed_Reads_$subdir/ -t 40 -type TruHT;\
   rm_shorter.pl -f 'Adpt_Trimmed_Reads_'$subdir'/2_Trim_adaptor_kmer/*.gz' -t 40;\
   rm Adpt_Trimmed_Reads_$subdir/1_Trim_garole/ -r;\
done;

