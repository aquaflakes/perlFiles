
#!/bin/bash

for i in 0deg 10deg 20deg 30deg 40deg 50deg H2O-25deg 40perEG-m20deg 40perD2O-25deg 80perD2O-25deg 40perEG-25deg 80perEG-25deg;
do	

ls ~/Nut_zhuData/seqFiles2/FJ10tempSolVar/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1III-FJ10tempSolvVar-${i}*IIIc4*.gz | \
parallel -j 20 gkmer_cnt.R -f {} -k 3 -o gkmer_cnt/gkmer_cnt_${i}c4 -t 20 --rmdup

done;
