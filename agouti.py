import os
import sys
import re
import argparse
import logging

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

from src import agouti_sequence as agSeq
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

	exclusiveGroup = parser.add_mutually_exclusive_group(required=True)
	exclusiveGroup.add_argument("-contig",
						metavar="FILE",
						dest="contigFasta",
						help="specify the assembly in contigs in FASTA format")
	exclusiveGroup.add_argument("-scaffold",
						metavar="FILE",
						dest="scaffoldFasta",
						help="specify the assembly in scaffolds in FASTA format")
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
	parser.add_argument("-outdir",
						metavar="DIR",
						dest="outDir",
						default=".",
						required=True,
						help="specify the directory to store output files")
	parser.add_argument("-p",
						metavar="STR",
						dest="prefix",
						default="agouti",
						help="specify the prefix for all output files [agouti]")
	parser.add_argument("-mnl",
						metavar="INT",
						dest="minSupport",
						default=5,
						help="minimum number of joining reads pair supports [5]")
	parser.add_argument("-nN",
						metavar="INT",
						dest="numNs",
						default=1000,
						help="number of Ns put in between a pair of contigs [1000]")
	parser.add_argument("-debug",
						action='store_true',
						help="specify the output prefix")

	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	bamFile = args.bamFile
	gffFile = args.gff
	prefix = args.prefix
	outDir = os.path.realpath(args.outDir)
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	if args.contigFasta:
		seqNames, dSeq = agSeq.get_contigs(args.contigFasta)
	elif args.scaffoldFasta:
		seqNames, dSeq = agSeq.get_scaffolds(args.scaffoldFasta)

	dGFFs = agGFF.get_gene_models(gffFile)

	dContigPairs = agBAM.get_joining_pairs(bamFile, args.minSupport)

	joinPairsFile = os.path.join(outDir, "%s.join_pairs" %(prefix))
#	dCtgPair2GenePair = agFILTER.map_contigPair2genePair(dContigPairs, dGFFs, joinPairsFile, args.minSupport)
	dCtgPair2GenePair = agFILTER.cleanContigPairs(dContigPairs, dGFFs, joinPairsFile, args.minSupport)

	pathList, edgeSenseDict, visitedDict = agSCAFF.agouti_scaffolding(seqNames, joinPairsFile, args.minSupport)
	agUPDATE.agouti_update(pathList, dSeq, seqNames,
						   edgeSenseDict, visitedDict, dGFFs,
						   dCtgPair2GenePair, outDir, prefix, args.numNs)

if __name__ == "__main__":
	main()
