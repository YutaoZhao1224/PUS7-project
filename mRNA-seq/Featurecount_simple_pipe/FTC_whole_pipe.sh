### usage: bash ~/Genome/tools_custom/RNA-seq/Featurecount_simple_pipe/FTC_whole_pipe.sh *shNC*bam *shFTO*bam shNC 2 shFTO 2
### while using paired-end read: add -p -B parameters

featureCounts -a ~/Genome/GTF/Gencode-V42-hg38-mRNA.gtf -s 2 \
	-T 20 -t exon -g gene_name -p -B --countReadPairs \
	-o mRNA-${@:$#-3:1}_vs_${@:$#-1:1}.txt ${@:1:$#-4} 2> FTC.log
grep 'Successfully assigned alignments' FTC.log | awk '{print $6}' > reads_num.txt
python ~/Genome/tools_custom/RNA-seq/Featurecount_simple_pipe/FeatureCounts_Asistant.py mRNA-${@:$#-3:1}_vs_${@:$#-1:1}.txt reads_num.txt ${@:$#-3:4}
rm FTC.log mRNA-${@:$#-3:1}_vs_${@:$#-1:1}.txt reads_num.txt mRNA-${@:$#-3:1}_vs_${@:$#-1:1}.txt.summary





#for i
#       	in "$@"; do
#	echo $i
#done

#echo ${@:$#-3:4}

#echo ${@:2:$#-5}
#echo mRNA-${@:$#-3:1}_vs_${@:$#-1:1}.txt
#echo $#

