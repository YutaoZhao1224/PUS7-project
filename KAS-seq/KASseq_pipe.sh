R1=$1
R2=$2

for i in ./log ./mapping ./mapping/bw
do
	if [ ! -d $i ];then
		mkdir $i
	fi
done

cutadapt -j 15 -m 26 --max-n=0 -e 0.15 -q 20 --nextseq-trim=20 -O 6 -n 3 --pair-filter=any \
	-a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGT \
	-o ${R1/.fastq.gz/.trimmed.fastq.gz} -p ${R2/.fastq.gz/.trimmed.fastq.gz} \
	$R1 $R2 > ./log/${R1/.R1.fastq.gz/}.log
clumpify.sh in=${R1/.fastq.gz/.trimmed.fastq.gz} in2=${R2/.fastq.gz/.trimmed.fastq.gz} \
	out=${R1/.fastq.gz/.dedupe.fastq.gz} out2=${R2/.fastq.gz/.dedupe.fastq.gz} dedupe \
	2>> ./log/${R1/.R1.fastq.gz/}.log

### remove redundant files
rm ${R1/.fastq.gz/.trimmed.fastq.gz} ${R2/.fastq.gz/.trimmed.fastq.gz}

### mapping to hg38 genome
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 0 -X 1000 --mp 5 --un-conc-gz \
	${R1/.R1.fastq.gz/}.fastq.gz -x ~/Genome/hg38_UCSC \
	-1 ${R1/.fastq.gz/.dedupe.fastq.gz} -2 ${R2/.fastq.gz/.dedupe.fastq.gz} \
	-S ./mapping/${R1/.R1.fastq.gz/}.mapped.sam 2>> ./log/${R1/.R1.fastq.gz/}.log
awk '$5 >= 20 || $1 ~ /^@/' ./mapping/${R1/.R1.fastq.gz/}.mapped.sam | samtools sort -@ 10 \
	--input-fmt-option "filter=(flag == 99 || flag == 147 || flag == 83 || flag == 163 )" \
	-O BAM -o ./mapping/${R1/.R1.fastq.gz/}.sorted.bam
samtools index -@ 6 ./mapping/${R1/.R1.fastq.gz/}.sorted.bam
rm ./mapping/${R1/.R1.fastq.gz/}.mapped.sam ${R1/.fastq.gz/.dedupe.fastq.gz} ${R2/.fastq.gz/.dedupe.fastq.gz} \
	${R1/.R1.fastq.gz/}.fastq.1.gz ${R1/.R1.fastq.gz/}.fastq.2.gz

### make bigwig file
bamCoverage -b ./mapping/${R1/.R1.fastq.gz/}.sorted.bam --samFlagInclude 128 --outFileFormat bigwig \
	-p 10 -bl ~/Genome/blacklist/hg38-blacklist.v2.bed --effectiveGenomeSize 3113504866 \
	--normalizeUsing RPKM -o ./mapping/bw/${R1/.R1.fastq.gz/}.bw


