import os
import sys
import re
import collections
import time

from lib import agouti_log as agLOG

def getCIGAR(cigar):
	tmp_cigar = re.split("([MIDNSHPX=])", cigar)[:-1]
	alnLen = 0.0
	for i in range(0, len(tmp_cigar), 2):
		if tmp_cigar[i+1] == 'M':
			alnLen += int(tmp_cigar[i])
	return alnLen

def explainSAMFlag(samFlag):
	bits = []
	for i in range(12):
		bits.append(2**(12-i-1))
	paired, proper = 0, 0
	selfUnMapped, mateUnMapped = 0, 0
	selfStrand, mateStrand = "+", "+"
	secondAln = 0
	duplicates = 0
	for i in range(len(bits)):
		if samFlag >= bits[i]:
			if bits[i] == 1:
				paired = 1
			elif bits[i] == 2:
				proper = 1
			elif bits[i] == 4:
				selfMapped = 0
			elif bits[i] == 8:
				mateMapped = 1
			elif bits[i] == 16:
				selfStrand = "-"
			elif bits[i] == 32:
				mateStrand = "-"
			elif bits[i] == 256:
				secondAln = 1
			elif bits[i] == 1024:
				duplicates = 1
			samFlag -= bits[i]

	return (paired, proper, selfUnMapped, mateUnMapped,
			selfStrand, mateStrand, secondAln, duplicates)

def getMismatches(tags):
	nMismatches = -1
	for i in range(len(tags)):
		tmp_tag = tags[i].split(':')
		if tmp_tag[0] == "NM":
			nMismatches = int(tmp_tag[2])
	return nMismatches

def getMappedRegionOnContigs(start, alnLen, flags):
	if flags[4] == "+":
		return start, start + alnLen
	else:
		return start+alnLen, start

def set_module_name(name):
	global moduleName
	moduleName = name

def try_continue_last_run(moduleProgressLogFile, moduleOutputFile):
	if os.path.exists(moduleProgressLogFile):
		fLOG = open(moduleProgressLogFile, 'r')
		lastLine = fLOG.readlines()[-1].strip()
		fLOG.close()
		moduleProgressLogger = moduleProgressLogObj.create_logger(moduleProgressLogFile)
		moduleProgressLogger.info("[BEGIN] Trying to read joining pairs from last run")
		moduleProgressLogger.info("Found log file from last run")
		moduleProgressLogger.info("Checking EXIT STATUS of last run")
		if lastLine.split('-')[-1].strip() == "Succeeded":
			moduleProgressLogger.info("Last run was successful")
			moduleProgressLogger.info("Skipping re-reading BAM file")
			moduleProgressLogger.info("[BEGIN] Reading output from last run")
			nJoinPairs = 0
			dContigPairs = collections.defaultdict(list)
			with open(moduleOutputFile, 'r') as fJOINPAIR:
				for line in fJOINPAIR:
					tmpLine = line.strip().split("\t")
					readsID = tmpLine[0]
					contigA, startA, stopA, senseA = tmpLine[1:5]
					contigB, startB, stopB, senseB = tmpLine[5:]
					startA = int(startA)
					stopA = int(stopA)
					startB = int(startB)
					stopB = int(stopB)
					nJoinPairs += 1
					if contigA <= contigB:
						if (contigA, contigB) not in dContigPairs:
							dContigPairs[contigA, contigB] = [(startA, startB, stopA, stopB, senseA, senseB, readsID)]
						else:
							dContigPairs[contigA, contigB] += [(startA, startB, stopA, stopB, senseA, senseB, readsID)]
					else:
						if (contigB, contigA) not in dContigPairs:
							dContigPairs[contigB, contigA] = [(startB, startA, stopB, stopA, senseB, senseA, readsID)]
						else:
							dContigPairs[contigB, contigA] += [(startB, startA, stopB, stopA, senseB, senseA, readsID)]
			moduleProgressLogger.info("%d joining pairs parsed" %(nJoinPairs))
			moduleProgressLogger.info("%d contig pairs given by these joining pairs" %(len(dContigPairs)))
			moduleProgressLogger.info("Succeeded")
			return dContigPairs, moduleProgressLogger
		else:
			moduleProgressLogger.info("Last run was NOT successful")
			return None, moduleProgressLogger
	return None, None

def get_joining_pairs(bamStream, moduleOutDir, prefix,
					  logLevel, overwrite, minMapQ=5,
					  minFracOvl=0.0, maxFracMismatch=1):
	moduleProgressLogFile = os.path.join(moduleOutDir, "%s.agouti_join_pairs.progressMeter" %(prefix))
	moduleDebugLogFile = os.path.join(moduleOutDir, "%s.agouti_join_pairs.debug" %(prefix))
	moduleOutputFile = os.path.join(moduleOutDir, "%s.agouti.join_pairs.unfiltered.txt" %(prefix))
	global moduleProgressLogObj
	moduleProgressLogObj = agLOG.AGOUTI_LOG(moduleName)
	moduleProgressLogger = None
	if not overwrite:
		dContigPairs, moduleProgressLogger = try_continue_last_run(moduleProgressLogFile, moduleOutputFile)
		if dContigPairs is not None:
			return dContigPairs
	moduleProgressLogger = moduleProgressLogObj.create_logger(moduleProgressLogFile)

	moduleDEBUGLogger = agLOG.AGOUTI_DEBUG_LOG(moduleName+"_DEBUG").create_logger(moduleDebugLogFile)

	with open(moduleOutputFile, 'w') as fOUT:
		moduleProgressLogger.info("[BEGIN] Identifying joining pairs")
		moduleProgressLogger.info("# processed\t| Current Reads ID\t| Elapsed Time")
		moduleDEBUGLogger.debug("Reads_ID\tLocationA\tLocationB\tmapQA\tmapQB\tsenseA\tsenseB\treadLenA\treadLenB")
		startTime = time.time()
		dContigPairs = collections.defaultdict(list)
		nJoinPairs = 0
		nReadsPairs = 0
		while True:
			pairA = bamStream.readline().strip().split("\t")
			pairB = bamStream.readline().strip().split("\t")
			# reach the end of the file
			if len(pairA) == 1 or len(pairB) == 1:
				break
			contigA = pairA[2]
			contigB = pairB[2]
			nReadsPairs += 1
			if pairA[0] == pairB[0] and contigA != contigB:
				readsID = pairA[0]
				alnLenA = getCIGAR(pairA[5])
				alnLenB = getCIGAR(pairB[5])
				leftMostPosA = int(pairA[3])
				leftMostPosB = int(pairB[3])
				readLenA = len(pairA[9])
				readLenB = len(pairB[9])
				nMismatchesA = getMismatches(pairA[11:])
				nMismatchesB = getMismatches(pairB[11:])
				mapQA = int(pairA[4])
				mapQB = int(pairB[4])
				flagsA = explainSAMFlag(int(pairA[1]))
				flagsB = explainSAMFlag(int(pairB[1]))
				senseA = flagsA[4]
				senseB = flagsB[4]
				moduleDEBUGLogger.debug("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d" %(readsID,
										contigA+":"+str(leftMostPosA), contigB+":"+str(leftMostPosB),
										mapQA, mapQB, senseA, senseB, readLenA, readLenB))

				if (min(alnLenA/readLenA, alnLenB/readLenB) >= minFracOvl and				# minimum fraction of overlaps
					max(nMismatchesA/alnLenA, nMismatchesB/alnLenB) <= maxFracMismatch and	# maximum fraction of mismatches
					min(mapQA, mapQB) >= minMapQ):				# minimum mapping quality
					startA = leftMostPosA + 1
					stopA = startA + 1 + int(alnLenA)
					startB = leftMostPosB + 1
					stopB = startB + 1 + int(alnLenB)
					nJoinPairs += 1
					if contigA <= contigB:
						if (contigA, contigB) not in dContigPairs:
							dContigPairs[contigA, contigB] = [(startA, startB, stopA, stopB, senseA, senseB, readsID)]
						else:
							dContigPairs[contigA, contigB] += [(startA, startB, stopA, stopB, senseA, senseB, readsID)]
						fOUT.write("%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\n" %(readsID, contigA, startA,
																			stopA, senseA, contigB,
																			startB, stopB, senseB))
					else:
						if (contigB, contigA) not in dContigPairs:
							dContigPairs[contigB, contigA] = [(startB, startA, stopB, stopA, senseB, senseA, readsID)]
						else:
							dContigPairs[contigB, contigA] += [(startB, startA, stopB, stopA, senseB, senseA, readsID)]
						fOUT.write("%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\n" %(readsID, contigB, startB,
																	stopB, senseB, contigA,
																	startA, stopA, senseA))
			if nReadsPairs % 5000000 == 0:
				elapsedTime = float((time.time() - startTime)/60)
				moduleProgressLogger.info("%d processed\t| %s\t| %.2f m" %(nReadsPairs, readsID, elapsedTime))

	moduleProgressLogger.info("%d joining pairs parsed" %(nJoinPairs))
	moduleProgressLogger.info("%d contig pairs given by these joining pairs" %(len(dContigPairs)))
	moduleProgressLogger.info("Succeeded")
	return dContigPairs
