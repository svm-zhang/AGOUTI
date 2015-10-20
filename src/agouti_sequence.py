import os
import sys
import logging

from lib import agouti_log as agLOG

def set_module_name(name):
	global moduleName
	moduleName = name

def get_contigs(contigFasta, moduleOutDir, prefix, logLevel):
	moduleLogFile = os.path.join(moduleOutDir, "%s.agouti_seq.progressMeter" %(prefix))
	moduleProgressLogger = agLOG.AGOUTI_LOG(moduleName).create_logger(moduleLogFile)
	try:
		fCONTIG = open(contigFasta, 'r')
	except IOError:
		moduleProgressLogger.error("Error opening contig file: %s" %(contigFasta), exc_info=True)
		sys.exit()

	moduleProgressLogger.info("[BEGIN] Reading the initial assembly")
	seq = ""
	contigs = []
	seqLens = []
	contigDict = {}
	contigIndex = 0
	with open(contigFasta, 'r') as fCONTIG:
		for line in fCONTIG:
			if line.startswith('>'):
				if seq != "":
					contigDict[contigIndex] = seq
					moduleProgressLogger.debug("%s\t%d" %(contig, len(seq)))
					seqLens.append(len(seq))
					contigIndex += 1
					seq = ""
				contig = line.strip()[1:]
				contigs.append(contig)
			else:
				seq += line.strip()
	# read one last sequence
	moduleProgressLogger.debug("%s\t%d" %(contig, len(seq)))
	contigDict[contigIndex] = seq
	seqLens.append(len(seq))

	n50 = get_assembly_NXX(seqLens)

	moduleProgressLogger.info("%d sequences parsed" %(len(contigDict)))
	moduleProgressLogger.info("The give assembly N50: %d" %(n50))
	moduleProgressLogger.info("[DONE]")

	return contigs, contigDict

def get_scaffolds(scaffoldFasta):
	moduleProgressLogger = agLOG.AGOUTI_LOG(moduleName, logLevel, None).create_logger()
	moduleProgressLogger.info("Processing scaffolds ... ")
	moduleProgressLogger.info("TO BE CONTINUED\n")

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
