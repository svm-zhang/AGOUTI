import os
import sys
import re
import argparse
import logging
import collections

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

from lib import agouti_log as agLOG
from src import agouti_sequence as agSeq
from lib import agouti_sam as agBAM
from src import agouti_filter_v2 as agFILTER
#from src import test_filter as agFILTER
from lib import agouti_gff as agGFF
from src import agouti_update as agUPDATE
from src import rnapathSTAR as agSCAFF

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
						type=int,
						default=5,
						help="minimum number of joining reads pair supports [5]")
	parser.add_argument("-nN",
						metavar="INT",
						dest="numNs",
						type=int,
						default=1000,
						help="number of Ns put in between a pair of contigs [1000]")
	parser.add_argument("-minMapQ",
						metavar="INT",
						dest="minMapQ",
						type=int,
						default=5,
						help="minimum of mapping quality to use [5]")
	parser.add_argument("-minFracOvl",
						metavar="FLOAT",
						dest="minFracOvl",
						type=float,
						default=0.0,
						help="minimum alignmentLen/readLen [0.0]")
	parser.add_argument("-maxFracMismatch",
						metavar="FLOAT",
						dest="maxFracMismatch",
						type=float,
						default=1.0,
						help="maximum fraction of mismatch of a give alignment [1.0]")
	parser.add_argument("-debug",
						action='store_true',
						help="specify the output prefix")
	parser.add_argument("-overwrite",
						action='store_true',
						help="specify whether to overwrite all results from last run [False]")

	args = parser.parse_args()
	return args

def main():

	args = parse_args()

	bamFile = args.bamFile
	gffFile = os.path.realpath(args.gff)
	prefix = args.prefix
	outDir = os.path.realpath(args.outDir)
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	logLevel = logging.INFO
	if args.debug:
		logLevel = logging.DEBUG
	mainLogFile = os.path.join(outDir, "%s.main.log" %(prefix))
	logger = agLOG.AGOUTI_LOG("Main").create_logger(mainLogFile)
	logger.info("Fasta File: %s" %(os.path.realpath(args.contigFasta)))
	logger.info("GFF file: %s" %(gffFile))
	logger.info("Output directory: %s" %(outDir))
	logger.info("Output prefix: %s" %(prefix))
	logger.info("Minimum number of supports: %d" %(args.minSupport))
	logger.info("Length of gaps filled: %d" %(args.numNs))

	agSeq.set_module_name("AGOUTI_Seq")
	moduleOutDir = os.path.join(outDir, "agouti_seq")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	if args.contigFasta:
		seqNames, dSeq = agSeq.get_contigs(args.contigFasta, moduleOutDir, prefix, logLevel)
	elif args.scaffoldFasta:
		seqNames, dSeq = agSeq.get_scaffolds(args.scaffoldFasta, logLevel)

	agGFF.set_module_name("AGOUTI_GFF")
	moduleOutDir = os.path.join(outDir, "agouti_GFFs")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	dGFFs = agGFF.get_gene_models(gffFile, moduleOutDir, prefix, logLevel)

	agBAM.set_module_name("AGOUTI_JoinPair")
	moduleOutDir = os.path.join(outDir, "agouti_join_pairs")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	dContigPairs = agBAM.get_joining_pairs(bamFile, moduleOutDir, prefix,
										   logLevel, args.overwrite,
										   args.minMapQ, args.minFracOvl,
										   args.maxFracMismatch)

#	joinPairsFile = os.path.join(outDir, "%s.join_pairs" %(prefix))
#	dCtgPair2GenePair = agFILTER.map_contigPair2genePair(dContigPairs, dGFFs, joinPairsFile, args.minSupport)
	agFILTER.set_module_name("AGOUTI_FILTER")
	# this module uses the same output directory as last module
	dCtgPair2GenePair, joinPairsFile = agFILTER.enforce_filters(dContigPairs, dGFFs, moduleOutDir,
																prefix, args.minSupport)

	agSCAFF.set_module_name("rnapathSTAR")
	moduleOutDir = os.path.join(outDir, "rnapathSTAR")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	scafPaths, edgeSenseDict, visitedDict = agSCAFF.rnapathSTAR(seqNames, joinPairsFile, moduleOutDir, prefix, args.minSupport)

	agUPDATE.set_module_name("AGOUTI_UPDATE")
	moduleOutDir = os.path.join(outDir, "agouti_update")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	agUPDATE.agouti_update(scafPaths, dSeq, seqNames,
						   edgeSenseDict, visitedDict, dGFFs,
						   dCtgPair2GenePair, outDir, prefix,
						   moduleOutDir, args.numNs)

if __name__ == "__main__":
	main()
