import os
import sys
import argparse
import subprocess as sp
import shlex
import re
import resource

agoutiBase = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.insert(1, agoutiBase)

__version__ = sp.check_output(shlex.split("git describe --tag --always --dirty"),
										  cwd=os.path.dirname(os.path.realpath(__file__)))

from lib import agouti_log as agLOG
from lib import agouti_sam as agBAM
from lib import agouti_gff as agGFF
from src import agouti_sequence as agSeq
from src import agouti_denoise as agDENOISE
from src import agouti_update as agUPDATE
from src import agouti_scaffolding as agSCAFF
from src import agouti_shred as agSHRED
from src import agouti_path as agPATH

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
	scafParser.add_argument("-no_update_gff",
							action="store_true",
							help="specify whether to update gff. EXPERIMENTAL")
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
							help=("minimum ratio of alignment length to the read length. "
								  "Specify 0.0 to turn off. [0.0]"
								 ))
	scafParser.add_argument("-maxFracMM",
							metavar="FLOAT",
							dest="maxFracMM",
							type=float,
							default=1.0,
							help=("maximum fraction of mismatches per alignment allowed. "
								  "Specify 1.0 to turn off. [1.0]"
								 ))
	scafParser.add_argument("-debug",
							action="store_true",
							help="specify to have info for debugging")
	scafParser.add_argument("-overwrite",
							action="store_true",
							help="specify whether to overwrite all results from last run")
	scafParser.set_defaults(func=run_scaffolder)

	usage = "Shredding genome assembly into contigs at gaps of a minimum length"
	shredParser = subparsers.add_parser("shred",
										help=usage)
	shredParser.add_argument("-assembly",
							 metavar="FILE",
							 dest="assemblyFile",
							 required=True,
							 help="specify the assembly in FASTA format. REQUIRED")
	shredParser.add_argument("-gff",
							 metavar="FILE",
							 dest="gffFile",
							 help="specify the annotation to shred")
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
	shredParser.set_defaults(func=run_shredder)

	usage = "update AGOUTI to the latest stable version"
	updateParser = subparsers.add_parser("update",
										 formatter_class=argparse.RawTextHelpFormatter,
										 help="\t"+usage)
	updateParser.set_defaults(func=update_local)

	if(len(sys.argv) == 1):
		# exit when no command provided
		parser.print_help()
		sys.exit(1)

	return(parser.parse_args())

def run_shredder(args):
	assemblyFile = args.assemblyFile
	gffFile = args.gffFile
	prefix = args.prefix
	minGap = args.minGap
	minCtgLen = args.minCtgLen

	agSHRED.agouti_shred_main(assemblyFile, gffFile, prefix,
							  minGap, minCtgLen)

def run_scaffolder(args):
	bamFile = args.bamFile
	gffFile = os.path.realpath(args.gff)
	prefix = args.prefix
	outDir = os.path.realpath(args.outDir)
	if(not os.path.exists(outDir)):
		os.makedirs(outDir)

	paraLogFile = os.path.join(outDir, "%s.parameters.txt" %(args.prefix))
	para = agLOG.PROGRESS_METER(parse_args.__name__)
	para.add_file_handler(paraLogFile)
	para.logger.info("Assembly: %s" %(os.path.realpath(args.assemblyFile)))
	para.logger.info("Gene Model: %s" %(gffFile))
	if(args.oriScafPath):
		para.logger.info("Original scaffold path: %s" %(args.oriScafPath))
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

	dCtgPair2GenePair, dCtgPairDenoise = agDENOISE.denoise_joining_pairs(dContigPairs, dGFFs,
																	  vertex2Name, outDir,
																	  prefix, args.minSupport,
																	  args.debug)

	agoutiPaths, dSenses = agSCAFF.run_scaffolding(vertex2Name, dCtgPairDenoise,
													   dCtgPair2GenePair, outDir, prefix,
													   args.minSupport, args.debug)

	if(args.oriScafPath):
		agoutiPaths, dCtgPair2GenePair, dSenses = agPATH.agouti_path_main(agoutiPaths, dSenses,
																		  vertex2Name, dGFFs,
																		  dCtgPair2GenePair,
																		  args.oriScafPath, outDir,
																		  prefix)

	agUPDATE.agouti_update(agoutiPaths, dSeqs, vertex2Name,
						   dSenses, dGFFs, dCtgPair2GenePair,
						   outDir, prefix,
						   args.nFills, args.debug, args.no_update_gff)

	para.logger.info("Peak memory use: %.5f GB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024*1024)))

def update_local(args):
	'''
		update to latest version
	'''
	version = agLOG.PROGRESS_METER("UPDATE")
	repoDir = os.path.dirname(os.path.realpath(__file__))
	# first check git availability
	checkGitVersion = "git --version"
	p = sp.Popen(shlex.split(checkGitVersion), stdout=sp.PIPE, stderr=sp.PIPE)
	pout, perr = p.communicate()
	if(p.returncode):
		version.logger.info("Please check your PATH for git")
		version.logger.info("Update unsuccessful")
		sys.exit(1)
	# Then compare local with remote
	version.logger.info("Checking available updates of AGOUTI")
	checkLocal = "git log -n 1 --pretty=\"%%H\""
	localVersion = sp.check_output(shlex.split(checkLocal), cwd=repoDir).strip()
#	if remoteVersion != localVersion:
	gitCmd = "git ls-remote origin"
	heads = sp.check_output(shlex.split(gitCmd), cwd=repoDir).split("\n")
	tags = []
	dVersions = {}
	for line in heads:
		if(line):
			tmpLine = line.strip().split("\t")
			if(re.search("refs/tag", tmpLine[1])):
				if(re.search("\^\{\}$", tmpLine[1])):
					dVersions[tmpLine[1].strip("^{}")] = tmpLine[0]
					continue
				else:
					dVersions[tmpLine[1]] = ""
				tags.append(tmpLine[1])
	latesTag = sorted(tags)[-1]
	latestHash = dVersions[latesTag]
	if(latestHash != localVersion):
		gitCmd = "git fetch --all"
		p = sp.Popen(shlex.split(gitCmd), stdout=sp.PIPE, stderr=sp.PIPE, cwd=repoDir)
		pout, perr = p.communicate()
		if(p.returncode):
			version.logger.error("git fetch error: %s" %(perr))
			sys.exit(1)
		gitCmd = "git checkout -q %s -b %s" %(latesTag, latesTag.split("/")[-1])
		p = sp.Popen(shlex.split(gitCmd), stdout=sp.PIPE, stderr=sp.PIPE, cwd=repoDir)
		pout, perr = p.communicate()
		if(p.returncode):
			version.logger.error("git checkout error: %s" %(perr))
			sys.exit(1)
		version.logger.info("Update successful")
		sys.exit(0)
	version.logger.info("Current version is the LATEST. No need to update")

def main():
	args = parse_args()
	args.func(args)

if __name__ == "__main__":
    main()
