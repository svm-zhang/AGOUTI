import os
import sys
import argparse
import logging
import collections

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

__version__ = "v0.2.1"

from lib import agouti_log as agLOG
from src import agouti_sequence as agSeq
from lib import agouti_sam as agBAM
from src import agouti_denoise as agDENOISE
from lib import agouti_gff as agGFF
from src import agouti_update as agUPDATE
from src import agouti_scaffolding as agSCAFF

def parse_args():
	use_message = '''
	Welcome to AGOUTI!\n
	'''

	parser = argparse.ArgumentParser(description=use_message)

#	exclusiveGroup = parser.add_mutually_exclusive_group(required=True)
	parser.add_argument("-assembly",
						metavar="FILE",
						dest="assemblyFile",
						required=True,
						help="specify the assembly in FASTA format")
#	exclusiveGroup.add_argument("-scaffold",
#						metavar="FILE",
#						dest="scaffoldFasta",
#						help="specify the assembly in scaffolds in FASTA format")
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
	parser.add_argument("-algorithm",
						metavar="STR",
						dest="algorithm",
						default="gene",
						help="specify scaffolding algorith: gene model or weight priority [gene]")
	parser.add_argument("-outdir",
						metavar="DIR",
						dest="outDir",
						default=".",
						required=True,
						help="specify the base directory to store all output files")
	parser.add_argument("-p",
						metavar="STR",
						dest="prefix",
						default="agouti",
						help="specify the prefix for all output files [agouti]")
	parser.add_argument("-k",
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
	parser.add_argument("-minMQ",
						metavar="INT",
						dest="minMQ",
						type=int,
						default=5,
						help="minimum of mapping quality to use [5]")
	parser.add_argument("-minFracOvl",
						metavar="FLOAT",
						dest="minFracOvl",
						type=float,
						default=0.0,
						help="minimum percentage of alignment length: alnLen/readLen [0.0]")
	parser.add_argument("-maxFracMM",
						metavar="FLOAT",
						dest="maxFracMM",
						type=float,
						default=1.0,
						help="maximum fraction of mismatch of a give alignment [1.0]")
	parser.add_argument("-debug",
						action='store_true',
						help="specify to have info for debugging")
	parser.add_argument("-overwrite",
						action='store_true',
						help="specify whether to overwrite all results from last run")
	parser.add_argument("-v",
						"--version",
						action="version",
						version="AGOUTI {version}".format(version=__version__))

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

	paraLogFile = os.path.join(outDir, "%s.parameters.txt" %(args.prefix))
	para = agLOG.PROGRESS_METER(parse_args.__name__)
	para.add_file_handler(paraLogFile)
	para.logger.info("Assembly: %s" %(os.path.realpath(args.assemblyFile)))
	para.logger.info("Gene Model: %s" %(gffFile))
	para.logger.info("Output directory: %s" %(outDir))
	para.logger.info("Output prefix: %s" %(prefix))
	para.logger.info("Minimum number of supports: %d" %(args.minSupport))
	para.logger.info("Length of gaps filled: %d" %(args.numNs))

#	if args.assemblyFile:
	vertex2Name, dSeq = agSeq.read_assembly(args.assemblyFile, outDir,
											prefix, args.debug)
#	elif args.scaffoldFasta:
#		vertex2Name, dSeq = agSeq.get_scaffolds(args.scaffoldFasta, logLevel)

	dGFFs = agGFF.get_gene_models(gffFile, outDir, prefix, args.debug)

	dContigPairs = agBAM.get_joining_pairs(bamFile, outDir, prefix,
										   args.overwrite, args.minMQ,
										   args.minFracOvl, args.maxFracMM,
										   args.debug)

	dCtgPair2GenePair, joinPairsFile = agDENOISE.denoise_joining_pairs(dContigPairs, dGFFs,
																	  vertex2Name, outDir,
																	  prefix, args.minSupport,
																	  args.debug)

	scafPaths, edgeSenseDict = agSCAFF.run_scaffolding(args.algorithm, vertex2Name, joinPairsFile,
												   dCtgPair2GenePair, outDir, prefix,
												   args.minSupport, args.debug)

	agUPDATE.agouti_update(scafPaths, dSeq, vertex2Name,
						   edgeSenseDict, dGFFs,
						   dCtgPair2GenePair, outDir, prefix,
						   args.numNs, args.debug)

if __name__ == "__main__":
	main()
