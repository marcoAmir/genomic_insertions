#!/bin/bash -e

large_del_dir=/cluster/u/amirma/Phenomic_data/aquatic_deletions_2017/largeDel

if [ $(hostname) != 'dev.stanford.edu' ];then
	echo -e "\n\terror: ssh to dev to excute (script calls featureBits)\n\n"
	exit
fi

if [ "$#" -ne 4 ]; then
	echo -e "\nCreating genomic regions for independent lesions\n  Usage:\n  $0 -l <assembly list comma delimited> \
	-c <assembly list comma delimited>\n\t-l list of assemblies with lesions\n\t-c list of conserved species\n"
	exit 1
fi

flag1=$1
list1=$2
flag2=$3
list2=$4
intacts=""
lesions=""

# parse input 
if [ "$#" -eq 4 ]; then
	if [ $flag1 = "-l" ] && [ $flag2 = "-c" ]; then
		echo $list1 | sed -e "s/,/\n/g" > assembly_lesions
		lesions=`echo $list1 | sed -e "s/,//g"`
		echo $list2 | sed -e "s/,/\n/g" > assembly_intact
		intacts=`echo $list2 | sed -e "s/,//g"`
	elif [ $flag1 = "-c" ] && [ $flag2 = "-l" ]; then
		echo $list1 | sed -e "s/,/\n/g" > assembly_intact
		intacts=`echo $list1 | sed -e "s/,//g"`
		echo $list2 | sed -e "s/,/\n/g" > assembly_lesions
		lesions=`echo $list2 | sed -e "s/,//g"`
	else
		echo -e "\nCreating genomic regions for independent lesions\n  Usage:\n  $0 -l <assembly list comma delimited> -c <optional:assembly list comma delimited>\n\t-l list of assemblies with lesions\n\t-c list of conserved species [optional]"
	fi
fi
	
awk -v large_del_dir=${large_del_dir} -v lesions_bed=${lesions}".bed" \
 -F'\t' '{if(NR==1){printf "featureBits hg38 "}; printf large_del_dir"/"$0".gapsInChains.bed "} \
 END {print " -bed="lesions_bed}' assembly_lesions | bash

awk -v large_del_dir=${large_del_dir} -v lesions_bed=${lesions}".bed" -v out_bed="lesions_"${lesions}"_"${intacts}".bed"\
 -F'\t' '{if(NR==1){printf "featureBits hg38 "lesions_bed" "}; printf  large_del_dir"/"$0".intact.bed "} \
 END {print " -bed="out_bed}' assembly_intact | bash

rm -rf assembly_lesions assembly_intact
