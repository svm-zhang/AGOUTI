#!/usr/bin/python

import os
import sys
import re
import argparse
import subprocess

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

from lib import agouti_sam as agBAM
from src import agouti_filter_v2 as agFILTER
#from src import test_filter as agFILTER
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
						dest="minSupport",
						default=5,
						help="minimum number of joining reads pair supporting a link between a contig pair")
	parser.add_argument("-nN",
						metavar="INT",
						dest="numNs",
						default=1000,
						help="number of Ns put in between a pair of contigs")
	parser.add_argument("-p",
						metavar="STR",
						dest="prefix",
						default="agouti",
						help="specify the output prefix")

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

	dContigPairs = agBAM.get_joining_pairs(bamFile, args.minSupport)

	joinPairsFile = os.path.join(outDir, "%s.join_pairs" %(prefix))
#	dCtgPair2GenePair = agFILTER.map_contigPair2genePair(dContigPairs, dGFFs, joinPairsFile, args.minSupport)
	dCtgPair2GenePair = agFILTER.cleanContigPairs(dContigPairs, dGFFs, joinPairsFile, args.minSupport)

	pathList, edgeSenseDict, visitedDict = agSCAFF.agouti_scaffolding(len(contigDict), nameList, joinPairsFile, args.minSupport)
	agUPDATE.agouti_update(pathList, contigDict, nameList,
						   edgeSenseDict, visitedDict, dGFFs,
						   dCtgPair2GenePair, outDir, prefix, args.numNs)

if __name__ == "__main__":
	main()
