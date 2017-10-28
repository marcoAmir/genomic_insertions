#!/bin/bash -e

# this script gets a text file listing paths to bed files, and a minimal coverage threshold
# it concatenates these bed files and computes an overlapping count on a genome (currently hg38)
# it outputs only the fragments that are covered by at least the selected threshold

# NOTE THAT: script is defaulted to hg38 as the reference genome

if [[ "$#" -ne 2 ]]; then
        echo -e "\nOutputs genomic-regions covered by at least the user \
		specified amount of fragments in bed\n  Usage:\n\t$0 bed_file_list cov_thres\n"
        exit 1
fi

bed_file_list=$1 
cov_thres=$2

rm -rf all_fragments.*

while read curr_chr
do
	curr_chr=`echo -e ${curr_chr} | cut -d" " -f1`
	echo -e "\n\t...collecting all fragments from ${curr_chr}"
	while read curr_bed
	do
		awk -v chr=${curr_chr} -F'\t' '{if($1==chr) print $0}' ${curr_bed} >> all_fragments.bed
	done < ${bed_file_list}
done < chrom.sizes.valid

# computing the coverage
echo -e "\n\t...computing genome coverage with bedtools (genomecov)"
bedtools genomecov -bga -split -i all_fragments.bed -g chrom.sizes.valid > all_fragments.cov.bed



