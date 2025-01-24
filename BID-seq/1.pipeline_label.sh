#!/bin/bash
SN=$1
if [ ! -d ./log ];then
        mkdir ./log
fi

cutadapt -j 15 --max-n=0 -e 0.1 -q 20 -m 32 -n 2 -O 8 --nextseq-trim=20 \
    -a ATCACGAGATCGGAAGAGCA -g TTCTACAGTCCGACGATC -o ${SN}.trimmed_1.fastq.gz ${SN}.fastq.gz > ./log/${SN}.report
cutadapt -j 15 -m 32 -n 2 -a AGATCGGAAGAGCA -g AGTCCGACGATC -O 10 \
    -o ${SN}.trimmed_2.fastq.gz ${SN}.trimmed_1.fastq.gz >> ./log/${SN}.report
clumpify.sh in=${SN}.trimmed_2.fastq.gz out=${SN}.dedupe.fastq.gz dedupe 2>> ./log/${SN}.report
cutadapt -j 15 -u 5 -u -5 -m 22 --rename='{id}_{cut_prefix}{cut_suffix} {comment}' -o ${SN}.barcoded.fastq.gz ${SN}.dedupe.fastq.gz
rm ${SN}.trimmed_1.fastq.gz ${SN}.trimmed_2.fastq.gz ${SN}.dedupe.fastq.gz

ulimit -n 1000000
STAR --runThreadN 20 --genomeDir ~/Genome/hg38_UCSC --readFilesIn <(gunzip -c ${SN}.barcoded.fastq.gz) \
    --alignEndsType Local --outFilterMatchNmin 0 --outFilterMismatchNmax 5 \
    --outFilterMismatchNoverLmax 0.2 --scoreDelOpen -2 --scoreDelBase -2 \
    --outFilterMultimapNmax 140 --outSAMmultNmax -1 \
    --outReadsUnmapped Fastx --limitBAMsortRAM 8000000000 --outSAMtype BAM Unsorted \
    --outFileNamePrefix ./mapping/${SN}_
rm ${SN}.barcoded.fastq.gz 
samtools sort -@ 10 --input-fmt-option 'filter=(flag == 0 || flag == 16)' -o ./mapping/${SN}.sorted.bam ./mapping/${SN}_Aligned.out.bam
samtools index -@ 5 ./mapping/${SN}.sorted.bam
rm ./mapping/${SN}_Aligned.out.bam
cat ./mapping/${SN}_Log.final.out >> ./log/${SN}.report

if [ ! -d ./mapping/STAR_log_file ];then
        mkdir ./mapping/STAR_log_file
fi
mv ./mapping/${SN}_Log.final.out ./mapping/${SN}_Log.out ./mapping/${SN}_Log.progress.out ./mapping/${SN}_SJ.out.tab ./mapping/STAR_log_file

if [ ! -d ./mapping/unmapped ];then
        mkdir ./mapping/unmapped
fi
mv ./mapping/${SN}_Unmapped.out.mate1 ./mapping/unmapped

if [ ! -d ./mapping/dedup_bam ];then
        mkdir ./mapping/dedup_bam
fi

if [ ! -d ./Del_detection ];then
        mkdir ./Del_detection
fi

#umi_tools dedup -I ./mapping/${SN}.sorted.bam -S ./mapping/dedup_bam/${SN}.dedup.bam >> ./log/${SN}.report
#rm ./mapping/${SN}.sorted.bam ./mapping/${SN}.sorted.bam.bai
samtools mpileup -Q 13 --input-fmt-option 'filter=(flag == 0 || flag == 16)' -d 0 -f ~/Genome/hg38_UCSC.fa \
    ./mapping/${SN}.sorted.bam | python ~/Genome/tools_custom/BID-seq/BID-seq_mpileup_T_deletion_Treated.py > ./Del_detection/${SN}_Del.txt

