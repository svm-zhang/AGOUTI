import os
import sys
import itertools
import re
import time

from lib import agouti_log as agLOG

def agouti_seq_main(assemblyFile, outDir, prefix, debug=0):
	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_seq")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_seq.progressMeter" %(prefix))
	agSeqProgress = agLOG.PROGRESS_METER(moduleName)
	agSeqProgress.add_console_handler()
	agSeqProgress.add_file_handler(progressLogFile)
	#contigs, dSeqs = get_contigs(assemblyFile, agSeqProgress)
	agSeqProgress.logger.info("[BEGIN] Reading the initial assembly")
	dSeqs = {}
	dHeaders = {}
	contigs = []
	seqLens = []
	seqIndex = 0
	for header, seq in read_fasta(assemblyFile):
		# split header on any non-alphabetic character
		# use only the first of the return list
		#header = re.split("\W+", header)[0]
		if header not in dHeaders:
			contigs.append(header)
			dSeqs[seqIndex] = seq
			seqIndex += 1
			dHeaders[header] = 1
			seqLens.append(len(seq))
		else:
			agSeqProgress.logger.error("AGOUTI found DUPLICATED header: %s" %(header))
			agSeqProgress.logger.error("QUIT")
			sys.exit(1)

	n50 = get_assembly_NXX(seqLens)

	agSeqProgress.logger.info("%d sequences parsed" %(len(dSeqs)))
	agSeqProgress.logger.info("The given assembly N50: %d" %(n50))
	agSeqProgress.logger.info("[DONE]")
	return contigs, dSeqs

def read_fasta(assemblyFile):
	'''
		Thanks to brentp and code monk
		sharing this elegant way to read
		Fasta file
	'''
	with open(assemblyFile, 'r') as fASSEMBLY:
		seqIter = (k[1] for k in itertools.groupby(fASSEMBLY, lambda line: line[0]== ">"))
		for header in seqIter:
			#header = re.split("\W+", header.next()[1:].strip())[0]
			header = header.next()[1:].strip().split()[0]
			seq = "".join(s.strip() for s in seqIter.next())
			yield header, seq

def rc_seq(seq):
	"""
		return the reverse and complementary of
		the give sequence
	"""
	alphabetsInString = "ACGTNRMKWYSBDHVTGCANYKMWRSVHDBacgtnrmkwysbdhvtgcanykmwrsvhdb"
	alphabets = { alphabetsInString[i]:alphabetsInString[i+15] for i
				  in range(60)
				  if i<15 or 30<=i<=44 }
	return "".join([alphabets[base] for base in seq[::-1]])

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
