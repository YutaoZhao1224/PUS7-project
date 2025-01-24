#!/bin/bash
SN=$1

if [ ! -d ./log ];then
    mkdir ./log
fi

total_read_num=$(bc <<< "$(zcat ${SN}.fastq.gz | wc -l)/4")
cutadapt -j 15 --max-n=0 -e 0.1 -q 20 -m 32 -n 3 -O 8 --nextseq-trim=20 \
    -a ATCACGAGATCGGAAGAGCA -g TTCTACAGTCCGACGATC -o ${SN}.trimmed_1.fastq.gz ${SN}.fastq.gz \
    > ./log/${SN}.report

cutadapt -j 15 -m 32 -n 4 -a AGATCGGAAGAGCA -g AGTCCGACGATC -O 10 \
    -o ${SN}.trimmed_2.fastq.gz ${SN}.trimmed_1.fastq.gz >> ./log/${SN}.report
rm ${SN}.trimmed_1.fastq.gz
After_cutting_adaptor_read_num=$(bc <<< "$(zcat ${SN}.trimmed_2.fastq.gz | wc -l)/4")

clumpify.sh in=${SN}.trimmed_2.fastq.gz out=${SN}.dedupe.fastq.gz dedupe 2>> ./log/${SN}.report
rm ${SN}.trimmed_2.fastq.gz

cutadapt -j 15 -u 5 -u -5 -m 22 --rename='{id}_{cut_prefix}{cut_suffix} {comment}' \
    -o ${SN}.barcoded.fastq.gz ${SN}.dedupe.fastq.gz >> ./log/${SN}.report
rm ${SN}.dedupe.fastq.gz

zcat ${SN}.barcoded.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",(NR+3)/4"-1");} else if(NR%4==2) print;}' > ${SN}.fa
seqkit seq -M 120 ${SN}.fa > ${SN}.barcoded.fa ### remove reads with length more than 120
rm ${SN}.barcoded.fastq.gz ${SN}.fa
After_dedupe_read_num=$(bc <<< "$(cat ${SN}.barcoded.fa | wc -l)/2")

### remove remove reads mapped to ribo-RNA
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 1 --mp 5 --un ${SN}.ribo-removed.fa \
    -x ~/Genome/hg38-45S/hg38_45S -f ${SN}.barcoded.fa > /dev/null 2>> ./log/${SN}.report
rm ${SN}.barcoded.fa
After_ribo_removal_read_num=$(bc <<< "$(cat ${SN}.ribo-removed.fa | wc -l)/2")

### mapping
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-gz ${SN}.unmap.fa.gz \
    -x ~/Genome/hg38_UCSC -f ${SN}.ribo-removed.fa -S ${SN}.sam 2>> ./log/${SN}.report
samtools sort --input-fmt-option 'filter=(flag == 0 || flag == 16)' -@ 6 -o ${SN}.aligned.sam ${SN}.sam
samtools sort --input-fmt-option 'filter=(flag == 0 || flag == 16)' -@ 6 -o ${SN}.aligned.bam ${SN}.sam
samtools index ${SN}.aligned.bam
rm ${SN}.sam ${SN}.ribo-removed.fa
rm ${SN}.unmap.fa.gz
Mapped_read_num=$(grep -v '^@' ${SN}.aligned.sam | wc -l)

### Peaks Clustering
~/workspace/PAR_CLIP/PARpipe/scripts/editPARalyzerINIfile.pl ~/workspace/PAR_CLIP/PARpipe/scripts/Default_PARalyzer_Parameters.ini \
    ${SN}.aligned.sam ~/Genome/hg38_UCSC.2bit > ${SN}.PARalyzer_Parameters.ini
~/workspace/PAR_CLIP/PARpipe/scripts/PARalyzer 12G ${SN}.PARalyzer_Parameters.ini >> ./log/${SN}.report
rm ${SN}_PARalyzer_Utilized.sam ${SN}.groups ${SN}.distribution ${SN}.aligned.sam ${SN}.PARalyzer_Parameters.ini

awk -F ',' 'BEGIN {OFS="\t"} FNR>1 {print $1,$3,$4,$5"="$8"="$7"="$10"="$11"="$12,$11,$2}' ${SN}.clusters \
    | sort -k 1,1V -k 2,2n > ${SN}.raw.bed
rm ${SN}.clusters

### Statistics
echo '' >> ./log/${SN}.report
echo '@Total Read number: ' $total_read_num >> ./log/${SN}.report
echo '@Read number after adaptor cutting: ' $After_cutting_adaptor_read_num >> ./log/${SN}.report
echo '@Read number after deduplication: ' $After_dedupe_read_num >> ./log/${SN}.report
echo '@Read number after rRNA removal: ' $After_ribo_removal_read_num >> ./log/${SN}.report
echo '@Mapped read number: ' $Mapped_read_num >> ./log/${SN}.report
echo '@Too-short-read ratio: ' $(bc <<< "scale=2;100-100*$After_cutting_adaptor_read_num/$total_read_num")'%' >> ./log/${SN}.report
echo '@Read Duplication Ratio: ' $(bc <<< "scale=2;100-100*$After_dedupe_read_num/$After_cutting_adaptor_read_num")'%' >> ./log/${SN}.report
echo '@Ribosome RNA contamination ratio: ' $(bc <<< "scale=2;100-100*$After_ribo_removal_read_num/$After_dedupe_read_num")'%' >> ./log/${SN}.report
echo '@Total mapped ratio: ' $(bc <<< "scale=2;100*$Mapped_read_num/$After_dedupe_read_num")'%' >> ./log/${SN}.report

echo ''
echo 'Total Read number: ' $total_read_num
echo 'Read number after adaptor cutting: ' $After_cutting_adaptor_read_num
echo 'Read number after deduplication: ' $After_dedupe_read_num
echo 'Read number after rRNA removal: ' $After_ribo_removal_read_num
echo 'Mapped read number: ' $Mapped_read_num
echo 'Too-short-read ratio: ' $(bc <<< "scale=2;100-100*$After_cutting_adaptor_read_num/$total_read_num")'%'
echo 'Read Duplication Ratio: ' $(bc <<< "scale=2;100-100*$After_dedupe_read_num/$After_cutting_adaptor_read_num")'%'
echo 'Ribosome RNA contamination ratio: ' $(bc <<< "scale=2;100-100*$After_ribo_removal_read_num/$After_dedupe_read_num")'%'
echo 'Total mapped ratio: ' $(bc <<< "scale=2;100*$Mapped_read_num/$After_dedupe_read_num")'%'


### Peaks Annotation
cat ${SN}.raw.bed | python ~/Genome/tools_custom/PAR_CLIP/Clip_Annotation.py ${SN}.unid.bed > ${SN}.annotated.txt



