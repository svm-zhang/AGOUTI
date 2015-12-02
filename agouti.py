import os
import sys
import argparse
import subprocess as sp
import shlex
import re

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

__version__ = sp.check_output(shlex.split("git -C %s describe --tag --always --dirty"
										  %(os.path.dirname(os.path.realpath(__file__)))))

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

	parser.add_argument("--version",
						action="version",
						version="AGOUTI {version}".format(version=__version__))

	subparsers = parser.add_subparsers(title="Commands",
									   metavar="",
									   dest="command")

	usage = "Scaffolding genome assembly and update genome annotation"
	scafParser = subparsers.add_parser("scaffold",
									   formatter_class=argparse.RawTextHelpFormatter,
									   description=usage,
									   help="\t"+usage)
	scafParser.add_argument("-assembly",
							metavar="FILE",
							dest="assemblyFile",
							required=True,
							help="specify the assembly in FASTA format")
	scafParser.add_argument("-bam",
							metavar="FILE",
							dest="bamFile",
							default="-",
							required=True,
							help="specify the RNA-seq mapping results in BAM format")
	scafParser.add_argument("-gff",
							metavar="FILE",
							dest="gff",
							required=True,
							help="specify the predicted gene model in GFF format")
	scafParser.add_argument("-shredpath",
							metavar="FILE",
							dest="oriScafPath",
							help=("specify the original scaffolding path "
								  "obtained from genome shredding"))
	scafParser.add_argument("-outdir",
							metavar="DIR",
							dest="outDir",
							default=".",
							required=True,
							help="specify the base directory to store all output files")
	scafParser.add_argument("-p",
							metavar="STR",
							dest="prefix",
							default="agouti",
							help="specify the prefix for all output files [agouti]")
	scafParser.add_argument("-k",
							metavar="INT",
							dest="minSupport",
							type=int,
							default=5,
							help="minimum number of joining reads pair supports [5]")
	scafParser.add_argument("-nf",
							metavar="INT",
							dest="nFills",
							type=int,
							default=1000,
							help="number of Ns put in between a pair of contigs [1000]")
	scafParser.add_argument("-minMQ",
							metavar="INT",
							dest="minMQ",
							type=int,
							default=5,
							help="minimum of mapping quality to use [5]")
	scafParser.add_argument("-minFracOvl",
							metavar="FLOAT",
							dest="minFracOvl",
							type=float,
							default=0.0,
							help="minimum percentage of alignment length: alnLen/readLen [0.0]")
	scafParser.add_argument("-maxFracMM",
							metavar="FLOAT",
							dest="maxFracMM",
							type=float,
							default=1.0,
							help="maximum fraction of mismatch of a give alignment [1.0]")
	scafParser.add_argument("-debug",
							action="store_true",
							help="specify to have info for debugging")
	scafParser.add_argument("-overwrite",
							action="store_true",
							help="specify whether to overwrite all results from last run")
	scafParser.add_argument("--justrun",
						action="store_true",
						help="specify to just run AGOUTI without updating to latest version")
	scafParser.set_defaults(func=run_scaffolder)

	usage = "Shredding genome assembly into contigs at gaps of a minimum length"
	shredParser = subparsers.add_parser("shred",
										help=usage)
	shredParser.add_argument("-assembly",
							 metavar="FILE",
							 dest="assemblyFile",
							 required=True,
							 help="specify the assembly in FASTA format. REQUIRED")
	shredParser.add_argument("-p",
							 metavar="STR",
							 dest="prefix",
							 required=True,
							 help="specify a output prefix. REQUIRED")
	shredParser.add_argument("-mlg",
							 metavar="INT",
							 dest="minGap",
							 type=int,
							 default=5,
							 help="specify a minimum length of gaps to shred [5]")
	shredParser.add_argument("-mlc",
							 metavar="INT",
							 dest="minCtgLen",
							 type=int,
							 default=1000,
							 help="specify a minimum length of contigs [1000]")
	shredParser.add_argument("--justrun",
						action="store_true",
						help="specify to just run AGOUTI without updating to latest version")
	shredParser.set_defaults(func=run_shredder)

	if len(sys.argv) == 1:
		# exit when no command provided
		parser.print_help()
		sys.exit(1)

	return parser.parse_args()

def run_shredder(args):
	assemblyFile = args.assemblyFile
	prefix = args.prefix
	minGap = args.minGap
	minCtgLen = args.minCtgLen

	agSeq.assembly_breaker(assemblyFile, prefix,
						   minGap, minCtgLen)

def run_scaffolder(args):
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
	para.logger.info("Minimum number of supports: %d"
					 %(args.minSupport))
	para.logger.info("Length of gaps to fill between contigs: %d"
					 %(args.nFills))

	vertex2Name, dSeqs = agSeq.agouti_seq_main(args.assemblyFile, outDir,
											   prefix, args.debug)

	dGFFs = agGFF.get_gene_models(gffFile, outDir, prefix, args.debug)

	dContigPairs = agBAM.agouti_sam_main(bamFile, outDir, prefix,
										   args.overwrite, args.minMQ,
										   args.minFracOvl, args.maxFracMM,
										   args.debug)

	dCtgPair2GenePair, joinPairsFile = agDENOISE.denoise_joining_pairs(dContigPairs, dGFFs,
																	  vertex2Name, outDir,
																	  prefix, args.minSupport,
																	  args.debug)

	scafPaths, edgeSenseDict = agSCAFF.run_scaffolding(vertex2Name, joinPairsFile,
													   dCtgPair2GenePair, outDir, prefix,
													   args.minSupport, args.debug)

	agUPDATE.agouti_update(scafPaths, dSeqs, vertex2Name,
						   edgeSenseDict, dGFFs,
						   dCtgPair2GenePair, outDir, prefix,
						   args.oriScafPath,
						   args.nFills, args.debug)

def check_version(repoDir):
	pass

def update_local(version, repoDir):
	'''
		update to latest version
	'''
	# first check git availability
	checkGitVersion = "git --version"
	p = sp.Popen(shlex.split(checkGitVersion), stdout=sp.PIPE, stderr=sp.PIPE)
	pout, perr = p.communicate()
	if p.returncode:
		version.logger.info("Please check your PATH for git")
		version.logger.info("Skip trying to update AGOUTI")
		return
	# Then compare local with remote
	version.logger.info("Checking available updates of AGOUTI")
	checkRemote = "git -C %s ls-remote origin master" %(repoDir)
	output = sp.check_output(shlex.split(checkRemote))
	remoteVersion = output.strip().split("\t")[0]
	checkLocal = "git -C %s log -n 1 --pretty=\"%%H\"" %(repoDir)
	localVersion = sp.check_output(shlex.split(checkLocal)).strip()
	if remoteVersion != localVersion:
		gitCmd = "git -C %s ls-remote origin" %(repoDir)
		heads = sp.check_output(shlex.split(gitCmd)).split("\n")
		tags = []
		for line in heads:
			if line:
				tmpLine = line.strip().split("\t")
				if re.search("refs/tag", tmpLine[1]):
					if re.search("\^\{\}$", tmpLine[1]):
						continue
					tags.append(tmpLine[1])
		latesTag = sorted(tags)[-1]
		gitCmd = "git -C %s fetch --all" %(repoDir)
		p = sp.Popen(shlex.split(gitCmd), stdout=sp.PIPE, stderr=sp.PIPE)
		pout, perr = p.communicate()
		if p.returncode:
			version.logger.error("git fetch error: %s" %(perr))
			sys.exit(1)
		gitCmd = "git -C %s checkout -q %s -b %s" %(repoDir, latesTag, latesTag.split("/")[-1])
		p = sp.Popen(shlex.split(gitCmd), stdout=sp.PIPE, stderr=sp.PIPE)
		pout, perr = p.communicate()
		if p.returncode:
			version.logger.error("git checkout error: %s" %(perr))
			sys.exit(1)
		version.logger.info("Update successful")
		version.logger.info("Rebootting AGOUTI")
		main(0)
		sys.exit(0)
	return
#		version.logger.info("Please re-run AGOUTI with the latest version")
#

def main(update):
	args = parse_args()
	version = agLOG.PROGRESS_METER("MAIN")
	version.logger.info("update: %d" %(update))
	if update and not args.justrun:
		repoDir = os.path.dirname(os.path.realpath(__file__))
		update_local(version, repoDir)
#		if check_version(repoDir):
#	if check_version():
#		if not args.justrun:
#			update_local(version)
#			sys.exit(0)
	version.logger.info("Running AGOUTI with current version: %s" %(__version__.strip()))
	args.func(args)

if __name__ == "__main__":
	main(1)
