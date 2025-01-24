#!/bin/bash

bam_file=$1
if [[ $bam_file == *input* ]]; then
	samtools mpileup -Q 20 --input-fmt-option 'filter=(flag == 0 || flag == 16)' -d 0 -f ~/Genome/hg38_UCSC.fa \
		./mapping/${SN}.sorted.bam | python ~/Genome/tools_custom/BID-seq/BID-seq_mpileup_T_deletion_Input.py \
		> ./Del_detection/${SN}_Del.txt
fi

if [[ $bam_file == *label* ]]; then
	samtools mpileup -Q 20 --input-fmt-option 'filter=(flag == 0 || flag == 16)' -d 0 -f ~/Genome/hg38_UCSC.fa \
		./mapping/${SN}.sorted.bam | python ~/Genome/tools_custom/BID-seq/BID-seq_mpileup_T_deletion_Input.py \
		> ./Del_detection/${SN}_Del.txt
fi


