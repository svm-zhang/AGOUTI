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

def agouti_shred_main(assemblyFile, gffFile, prefix,
					  minGaps, minCtgLen):
	pass

def assembly_breaker(assemblyFile, prefix, minGaps, minCtgLen):
	breakerProgress = agLOG.PROGRESS_METER("SHREDDER")
	breakerProgress.logger.info("[BEGIN] Shredding assembly")

	outdir = os.path.dirname(os.path.realpath(prefix))
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	outDebugFile = prefix + ".shred.debug"
	breakDebug = agLOG.DEBUG("SHREDDER", outDebugFile)
	outFa = prefix + ".ctg.fasta"
	outInfo = prefix +".shred.info.txt"
	with open(outFa, 'w') as fOUTFA, open(outInfo, 'w') as fINFO:
		genomeSize = 0
		splitSize = 0
		numContigs = 0
		contigLens = []
		nSeqs = 0
		startTime = time.time()
		breakerProgress.logger.info("# processed\t| Current sequence ID\t| Elapsed Time")
		for header, seq in read_fasta(assemblyFile):
			nSeqs += 1
			breakDebug.debugger.debug(">%s" %(header))
			genomeSize += len(seq)
			gapIndices = [(m.start(), m.end()) for m in re.finditer("[N|n]{%d,}" %(minGaps), seq)]
			gapIndices.append((len(seq), -1))
			breakDebug.debugger.debug("gapIndices: %s" %(str(gapIndices)))
			gapLens = []
			intervals = []
			if len(gapIndices) == 1:
				intervals.append((0, gapIndices[0][0]))
			elif gapIndices[-1][0] < minCtgLen:
				intervals.append((0, gapIndices[-1][0]))
			else:
				start = 0
				i = 0
				for i in range(len(gapIndices)):
					stop = gapIndices[i][0]
					breakDebug.debugger.debug("start %d stop %d" %(start, stop))
					if gapIndices[len(gapIndices)-1][0]-start < minCtgLen and \
					   len(gapIndices) > 1:
						breakDebug.debugger.debug("last short")
						breakDebug.debugger.debug("gapIndices[i]: %s" %(str(gapIndices[i])))
						breakDebug.debugger.debug("intervals: %s" %(intervals))
						#if len(intervals) > 0:
						intervals[-1] = (intervals[-1][0], gapIndices[len(gapIndices)-1][0])
						#else:
						#	intervals.append((start, gapIndices[len(gapIndices)-1][0]))
						break
					if stop-start+1 < minCtgLen:
						breakDebug.debugger.debug("short")
						breakDebug.debugger.debug("previous %s next %s" %(str(gapIndices[i-1]), str(gapIndices[i])))
						breakDebug.debugger.debug("length: %d" %(stop-start+1))
						continue
					if i < len(gapIndices)-1:
						gapLens.append(gapIndices[i][1]-gapIndices[i][0]+1)
					intervals.append((start, stop))
					start = gapIndices[i][1]+1
			breakDebug.debugger.debug("intervals: %s" %(intervals))
			breakDebug.debugger.debug("gapLen: %s" %(str(gapLens)))
			contigs = []
			for i in range(len(intervals)):
				start = intervals[i][0]
				stop = intervals[i][1]
				splitSize += (stop-start)
				contigID = "%s_%d" %(header, i)
				contigs.append(contigID)
				contigLens.append(stop-start)
				fOUTFA.write(">%s\n%s\n" %(contigID, seq[start:stop]))
			numContigs += len(contigs)
			if nSeqs % 10000 == 0:
				elapsedTime = float((time.time() - startTime)/60)
				breakerProgress.logger.info("%d processed\t| %s\t | %.2f m" %(nSeqs, header, elapsedTime))
			fINFO.write(">%s\n" %(header))
			if len(contigs) == 1:
				fINFO.write("%s\tNA\tNA\n" %(contigs[0]))
				continue
			for i in range(1, len(contigs)):
				fINFO.write("%s\t%s\t%d\n" %(contigs[i-1], contigs[i], gapLens[i-1]))
		if nSeqs < 10000:
			elapsedTime = float((time.time() - startTime)/60)
			breakerProgress.logger.info("%d processed\t| %s\t | %.2f m" %(nSeqs, header, elapsedTime))
		n50 = get_assembly_NXX(contigLens)
		breakerProgress.logger.info("Total length of the given assembly: %d"
									%(genomeSize))
		breakerProgress.logger.info("Total length of the shred assembly: %d"
									%(splitSize))
		breakerProgress.logger.info("Number of sequences in the shred assembly: %d"
									%(numContigs))
		breakerProgress.logger.info("N50 of the shred assembly: %d" %(n50))
