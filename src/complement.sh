#!/bin/bash -e 

for curr_assembly in bosTau8 loxAfr3 mm10 galGal5 canFam3
do
	echo ${curr_assembly}
	featureBits hg38 ${curr_assembly}.gapsInChains.allIncluded.bed -bed=${curr_assembly}.intact.bed
done
