#!/usr/bin/env python
#
# Amir Marcovitz
# Date created: 10/26/2017

# This code scans the chain files for a selected assembly and mark all the alignable regions (in output bed file)
# which are simply the coordinates of gapless blocks in the chains

import re
import sys
import argparse
import os
import gzip
import subprocess

def readArgs():
	print ""
	parser = argparse.ArgumentParser('Generate bed files and browser visualization for mammalian insertions')
	parser.add_argument('assembly', help='target assembly')
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	return args

def readAlignedRegions(assembly, subsetChains=True, refAssembly="hg38"):
	print("\n\t...extracting chain blocks for target assembly " + assembly + "\n")
	chainMapDir = "/cluster/u/amirma/geneLoss/" + refAssembly + "/mappings/"
	chainFile='/cluster/u/amirma/Phenomic_data/aquatic_deletions_2017/chain_files/' + refAssembly + '.' + assembly + '.all.chain.gz'
	outBed = open('./alignables_bed/' + assembly + '.blocksInChains.bed', 'w')
	refChr = 'chr1' 								# just to pass the first control
	chainIDs = GetSubsetChains(assembly, chainMapDir, chainFile, subsetChains)	# an array with chain IDs (human to query to analyze) 
	l=1
	with gzip.open(chainFile) as f:
		for line in f:
			print l
			l += 1
			if ((len(line.split('\t'))<=1 and len(line.split(' '))<=1) or line[0]=='#'):
				continue
			if re.search('chain', line):			# new chain
				refChr  = line.split(' ')[2]		# chr in reference genome
				chainID = line.split(' ')[12].rstrip()	# the id of the chain
				if (not conventionalChromosomeHG(refChr)) or (chainID not in chainIDs):
					while line != "\n":
						line = next(f)
						l += 1
					continue
				refPos  = int(line.split(' ')[5])	# the start pos of the chain - on reference
				refPos_n = refPos
				blockInd  = 0				# index of the block
			else:
				blockInd += 1
				refPos_n = refPos_n + int(line.split('\t')[0])		# Propagate the ref position to the end of current alignment block
				if int(line.split('\t')[1])>0 or len(line.split('\t'))==1:
					outBed.write(refChr + '\t' + str(refPos) + '\t' + str(refPos_n) + '\tblock' + str(blockInd) + '_' + str(chainID) + '\n')
					refPos = refPos_n
					refDelSize = int(line.split('\t')[1])
					refPos  = refPos + refDelSize
					refPos_n = refPos
	outBed.close()

def GetSubsetChains(assembly, chainMapDir, chainFile, subsetChains=True):
	print assembly, subsetChains
	chainIDs = []
	if (not subsetChains):
		with gzip.open(chainFile) as f:
			for line in f:
				if re.search('chain', line):
					chainIDs.append(line.split(' ')[12].rstrip())
	else:
		with open(chainMapDir + assembly + '.chain_ids') as f:
			for line in f:
				chainIDs.append(line.split(' ')[1])
	return list(set(chainIDs))
	

def conventionalChromosomeHG(refChr):
	a = (re.search('chr\d+\Z', refChr) and int(re.search('chr(\d+)\Z', refChr).group(1))<=22 and int(re.search('chr(\d+)\Z', refChr).group(1))>=1)
	b = (refChr=='chrM') or (refChr=='chrX') or (refChr=='chrY')
	if (a or b):
		return 1
	else:
		return 0

if __name__ == "__main__":
	print len(sys.argv)
	args = readArgs()
	readAlignedRegions(**args)
