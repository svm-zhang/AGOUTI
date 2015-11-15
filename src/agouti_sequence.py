import os
import sys

from lib import agouti_log as agLOG

def get_contigs(assemblyFile, agSeqProgress):
#	try:
#		fCONTIG = open(assemblyFile, 'r')
#	except IOError:
#		agSeqProgress.logger.error("Error opening contig file: %s" %(assemblyFile), exc_info=True)
#		sys.exit()

	agSeqProgress.logger.info("[BEGIN] Reading the initial assembly")
	seq = ""
	contigs = []
	seqLens = []
	contigDict = {}
	contigIndex = 0
	with open(assemblyFile, 'r') as fCONTIG:
		for line in fCONTIG:
			if line.startswith('>'):
				if seq != "":
					contigDict[contigIndex] = seq
					agSeqProgress.logger.debug("%s\t%d" %(contig, len(seq)))
					seqLens.append(len(seq))
					contigIndex += 1
					seq = ""
				contig = line.strip()[1:]
				contigs.append(contig)
			else:
				seq += line.strip()
	# read one last sequence
	agSeqProgress.logger.debug("%s\t%d" %(contig, len(seq)))
	contigDict[contigIndex] = seq
	seqLens.append(len(seq))

	n50 = get_assembly_NXX(seqLens)

	agSeqProgress.logger.info("%d sequences parsed" %(len(contigDict)))
	agSeqProgress.logger.info("The given assembly N50: %d" %(n50))
	agSeqProgress.logger.info("[DONE]")

	return contigs, contigDict

def read_assembly(assemblyFile, outDir, prefix, debug=0):
	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_seq")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_seq.progressMeter" %(prefix))
	agSeqProgress = agLOG.PROGRESS_METER(moduleName)
	agSeqProgress.add_file_handler(progressLogFile)
	contigs, contigDict = get_contigs(assemblyFile, agSeqProgress)

	return contigs, contigDict

def get_assembly_NXX(seqLens, nXX=50):
	seqLenSum = sum(seqLens)

	nXXThreshold = seqLenSum * (nXX/100.0)

	seqLens.sort()
	cumuSeqLen = 0
	nXXLen = 0
	for i in range(len(seqLens)-1, -1, -1):
		cumuSeqLen += seqLens[i]
		if cumuSeqLen > nXXThreshold:
			nXXLen = seqLens[i]
			break
	return nXXLen
