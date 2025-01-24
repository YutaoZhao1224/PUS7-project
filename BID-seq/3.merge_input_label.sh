#!/bin/bash

### usage: bash 3.merge_input_label.sh caRNA-shctr-BID-seq.input_Del.txt caRNA-shctr-BID-seq.label_Del.txt
if [ ! -d ./Merged ];then
	mkdir ./Merged
fi

if [[ $1 == *input* ]];then
	    input=$1
fi
if [[ $1 == *label* ]];then
	    label=$1
fi
if [[ $2 == *input* ]];then
	    input=$2
fi
if [[ $2 == *label* ]];then
	    label=$2
fi

awk 'BEGIN {OFS="\t"} NR==FNR{a[$1,$2,$3]=$4"\t"$5"\t"$6;next} \
	($1,$2,$3) in a{print $1,$2,$3,$4,$5,$6,a[$1,$2,$3]}' \
	$label $input \
	| awk 'BEGIN {OFS="\t"} {print $1,$2-1,$2,$5":"$6":"$8":"$9,$7-$4,$3}' \
	> ./Merged/$(basename $label .label_Del.txt).pU_site.txt

