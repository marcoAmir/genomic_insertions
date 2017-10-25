#!/usr/bin/env python

# Amir Marcovitz
# Date created: 10/24/2017

# collecting all insertions from human-mammalians liftOver chains
# this is just a minor twick to the script: /cluster/u/amirma/geneLoss/hg38/browser_visuals/assemblyLoFtracks.py
# this script iterates over our list of transcripts and for a given target assembly generates bed files for all collected insertions

import sys
import os
import argparse
import re
import gzip
import subprocess


import readInsertionsFromChains as RINS

chainMapDir = "/cluster/u/amirma/geneLoss/hg38/mappings/"

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Generate bed files and browser visualization for mammalian insertions')
	parser.add_argument('assembly', help='target assembly')
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args
	
def main(assembly, subsetChains=True, refAssembly="hg38"):
	includeAll=False	# select 'True' is you want gaps to be written down regardless of seq-gaps (good for outgroups)
	if not os.path.exists('./insertions_bed'):
    		os.makedirs('./insertions_bed')
	if not os.path.exists('./insertions_bed/' + assembly + '.insertionsInChains.bed') or not os.path.exists('./insertions_bed/' + assembly + '.overlappingGaps'):
		if assembly=='galGal5' or assembly=='phyCat1' or assembly=='balAcu1' or assembly=='lipVex1':
			subsetChains=False
		if assembly=='phyCat1' or assembly=='balAcu1' or assembly=='lipVex1':
			refAssembly="hg19"
		RINS.readGapsInChainFile(assembly, subsetChains, includeAll, refAssembly)

if __name__ == "__main__":
	print len(sys.argv)
	args = readArgs()
	main(**args)
