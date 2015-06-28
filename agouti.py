#!/usr/bin/python

import os
import sys
import re
import argparse
import subprocess

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

from lib import agouti_sam as agBAM
from src import agouti_filter as agFILTER
from lib import agouti_gff as agGFF
from src import agouti_update as agUPDATE
from src import agouti_scaffolding as agSCAFF

def parse_args():
	use_message = '''
	Welcome to AGOUTI!\n
	'''

	parser = argparse.ArgumentParser(description=use_message)

	parser.add_argument("-contig",
						metavar="FILE",
						dest="contigFile",
						required=True,
						help="specify the initial assembly in FASTA format")
	parser.add_argument("-bam",
						metavar="FILE",
						type=argparse.FileType('r'),
						dest="bamFile",
						default="-",
						required=True,
						help="specify the RNA-seq mapping results in BAM format")
	parser.add_argument("-gff",
						metavar="FILE",
						dest="gff",
						required=True,
						help="specify the predicted gene model in GFF format")
	parser.add_argument("-out",
						metavar="DIR",
						dest="outDir",
						default=".",
						required=True,
						help="specify the directory to store output files")
	parser.add_argument("-mnl",
						metavar="INT",
						dest="min_nLinks",
						default=5,
						help="minimum number of reads supporting a link between a contig pair")
	parser.add_argument("-p",
						metavar="STR",
						dest="prefix",
						default="agouti",
						help="specify the output prefix")
#	parser.add_argument("-sampleID",
#						metavar="STR",
#						dest="sampleID",
#						default="agouti",
#						help="specify the sampleID as the prefix of output files")
#	parser.add_argument("-isSE",
#						dest="isSE",
#						action="store_true",
#						default=False,
#						help="turn on this option if input is single-end RNAseq reads")
#	parser.add_argument("-sspacePATH",
#						dest="sspacePATH",
#						default="../SSPACE-STANDARD-3.0_linux-x86_64",
#						help="specify the path to SSPACE")
#	parser.add_argument("-agoutiPATH",
#						dest="agoutiPATH",
#						default=".",
#						help="specify the path to AGOUTI")
	# other possible arguments include "-compressFQ","-doFASTQC","doRSeqQC"

	args = parser.parse_args()
	return args

def get_initial_assembly(contigFile):
	try:
		fASSEMBLY = open(contigFile)
	except IOError:
		print "Error opening contig file: %s" % contigFile

	sys.stderr.write("Reading the initial assembly ... ")
	seq = ""
	nameList = []
	origSize = []
	contigDict = {}
	contigNum = 0
	with open(contigFile, 'r') as fASSEMBLY:
		for line in fASSEMBLY:
			if line.startswith('>'):
				if seq != "":
					contigDict[contigNum] = seq
					origSize.append(len(seq))
					contigNum += 1
					seq = ""
				chrom = line.strip()[1:]
				nameList.append(chrom)
			else:
				seq += line.strip()

	contigDict[contigNum] = seq
	origSize.append(len(seq))
	fASSEMBLY.close()
	sys.stderr.write("%d sequences parsed\n" %(len(contigDict)))
	return nameList, contigDict, origSize

def main():
	args = parse_args()

	contigFile = args.contigFile
	bamFile = args.bamFile
	gffFile = args.gff
	prefix = args.prefix
	outDir = os.path.realpath(args.outDir)
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	nameList, contigDict, origSize = get_initial_assembly(contigFile)

	dGFFs = agGFF.get_gene_models(gffFile)

	dContigPairs = agBAM.get_joining_pairs(bamFile, args.min_nLinks)

	joinPairsFile = os.path.join(outDir, "%s.join_pairs" %(prefix))
	dCtgPair2GenePair = agFILTER.cleanContigPairs(dContigPairs, dGFFs, joinPairsFile)

	pathList, edgeSenseDict, visitedDict = agSCAFF.agouti_scaffolding(len(contigDict), nameList, joinPairsFile)
	agUPDATE.agouti_update(pathList, contigDict, nameList, origSize,
						   edgeSenseDict, visitedDict, dGFFs,
						   dCtgPair2GenePair, outDir, prefix)

if __name__ == "__main__":
	main()
