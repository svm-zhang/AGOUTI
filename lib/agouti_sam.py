import os
import sys
import re
import collections
import time
import subprocess as sp
import shlex
import multiprocessing as mp

from lib import agouti_log as agLOG

def check_jp_entirity(dContigPairs, moduleOutDir):
	with open(os.path.join(moduleOutDir, "haha.jp"), 'w') as fOUT:
		for k, v in dContigPairs.iteritems():
			for shit in v:
				if k[0] <= k[1]:
					fOUT.write("%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\n" %(shit[-1], k[0], shit[0],
																	shit[2], shit[4], k[1],
																	shit[1], shit[3], shit[5]))
				else:
					fOUT.write("%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\n" %(shit[-1], k[1], shit[1],
																	shit[3], shit[5], k[0],
																	shit[0], shit[2], shit[4]))

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

def retrieve_joininng_pairs(agBAMProgress, agBAMOutAllJoinPairs):
	#progressLogFile = agBAMProgress.logFile
	#with open(progressLogFile, 'r') as fLOG:
	#	reports = fLOG.readlines()
	#if len(reports) == 0:
	#	return None
	#lastLine = reports[-1].strip()
	#exitStatus = lastLine.split('-')[-1].strip()
	#agBAMProgress.logger.info("---------------------------------")
	#agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")
	#agBAMProgress.logger.info("Found progress file from previous run")
	#if exitStatus == "Succeeded":
	#	agBAMProgress.logger.info("EXIT STATUS of last run: %s" %(exitStatus))
	#	agBAMProgress.logger.info("Skip re-reading BAM file")
	#	agBAMProgress.logger.info("Read joining-pairs from previous run")
	try:
		nJoinPairs = 0
		dContigPairs = collections.defaultdict(list)
		with open(agBAMOutAllJoinPairs, 'r') as fJOINPAIR:
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
		agBAMProgress.logger.info("%d joining pairs parsed" %(nJoinPairs))
		agBAMProgress.logger.info("%d contig pairs given by these joining pairs" %(len(dContigPairs)))
		return dContigPairs, nJoinPairs
	except KeyboardInterrupt:
		agBAMProgress.logger.info("Extract Joining-pairs INTERRUPTED by Keyboard")
		sys.exit(1)
	#else:
	#	return None

def check_samtools(agBAMProgress):
	'''
		check samtool availability
	'''
	samtoolsTest = "samtools"
	try:
		p = sp.Popen(shlex.split(samtoolsTest), stdout=sp.PIPE, stderr=sp.PIPE)
		pout, perr = p.communicate()
	except (OSError, ValueError) as e:
		agBAMProgress.logger.error("samtools is not available")
		agBAMProgress.logger.error("Please check if samtools is in your PATH")
		sys.exit(1)

def run_samtools(bamFile, agBAMProgress):
	'''
		run samtools to extract joining-pairs
	'''
	samtoolsCMD = "samtools view -h -F3328 %s" %(bamFile)
	# now do not need try
	try:
		p = sp.Popen(shlex.split(samtoolsCMD), stdout=sp.PIPE, stderr=sp.PIPE, bufsize=1)
		nRecords = 1
		for line in p.stdout:
			if line.startswith('@'):
				continue
			else:
				if nRecords%2 == 1:
					record = line.strip()
				else:
					record += "\n" + line.strip()
					yield record
				nRecords += 1
		perr = p.communicate()[1]
		# capture errors such as truncated BAM file
		if perr:
			#agBAMProgress.logger.error("SAMTOOLs issued ERROR when reading the BAM/SAM file")
			#for line in perr.strip().split("\n"):
			#	agBAMProgress.logger.error(line.strip())
			sys.exit(1)
	except (OSError, ValueError) as e:
		# in cases where invalid arguments called with samtools
		#agBAMProgress.logger.error("Error running SAMTOOLs")
		#agBAMProgress.logger.error("%s" %(e))
		sys.exit(1)

def worker(job):
	bamFile, config = job
	print bamFile, config, mp.current_process().name
	moduleName, moduleOutDir, prefix, overwrite, minMapQ, minFracOvl, maxFracMismatch, debug = config
	agBAMOutbase = os.path.basename(bamFile).strip(".bam")
	agBAMOutJoinPairs_file = os.path.join(moduleOutDir, agBAMOutbase+".jp")
	print agBAMOutJoinPairs_file

	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_jp.progressMeter" %(agBAMOutbase))
	agBAMProgress = agLOG.PROGRESS_METER(moduleName)
	#agBAMProgress.add_file_handler(progressLogFile)
	agBAMDone = os.path.join(moduleOutDir, "{}.done".format(agBAMOutbase))
	#if not os.path.exists(agBAMDone):
	#	agBAMProgress.add_file_handler(progressLogFile)
	#	agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")
	#else:
	if os.path.exists(agBAMDone):
		if not overwrite:
			agBAMProgress.add_file_handler(progressLogFile, 'a')
			agBAMProgress.logger.info("")
			agBAMProgress.logger.info("---------------------------------")
			agBAMProgress.logger.info("Found joining-pairs from previous run")
			agBAMProgress.logger.info("Reading joining-pairs from {}".format(agBAMOutJoinPairs_file))
			dContigPairs, nJoinPairs = retrieve_joininng_pairs(agBAMProgress, agBAMOutJoinPairs_file)
			if dContigPairs is not None:
				return (dContigPairs, nJoinPairs, bamFile)
			else:
				agBAMProgress.logger.info("Fail to pick up results from the previous run")
				agBAMProgress.logger.info("Re-processing the BAM file")
		else:
			agBAMProgress.add_file_handler(progressLogFile)
			agBAMProgress.logger.info("Found joining-pairs from previous run")
			agBAMProgress.logger.info("Overwritting")
			agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")

	agBAMDebug = None
	if debug:
		debugLogFile = os.path.join(moduleOutDir, "%s.agouti_jp.debug" %(agBAMOutbase))
		agBAMDebug = agLOG.DEBUG(moduleName, debugLogFile)
		print debugLogFile
	print progressLogFile

	try:
		with open(agBAMOutJoinPairs_file, 'w') as fOUT:
			agBAMProgress.logger.info("# processed\t| Current Reads ID\t| Elapsed Time")
			if debug:
				agBAMDebug.debugger.debug("Reads_ID\tLocationA\tLocationB\tmapQA\tmapQB\tsenseA\tsenseB\treadLenA\treadLenB")
			startTime = time.time()
			dContigPairs = collections.defaultdict(list)
			nJoinPairs = 0
			nReadsPairs = 0
			#for record in run_samtools(bamFile, agBAMProgress):
			for record in run_samtools(bamFile, None):
				tmpRecord = record.split("\n")
				pairA = tmpRecord[0].split("\t")
				pairB = tmpRecord[1].split("\t")
				readsID = pairA[0]
				contigA = pairA[2]
				contigB = pairB[2]
				mateCtgB = pairA[6]
				mateCtgA = pairB[6]
				nReadsPairs += 1
				# the first contidition makes sure
				# single end BAM are gonna have zero
				# joining-pairs extracted
				if contigA == "*" or contigB == "*":
					continue
				if pairA[0] == pairB[0] and contigA != contigB:
					alnLenA = getCIGAR(pairA[5])
					alnLenB = getCIGAR(pairB[5])
					leftMostPosA = int(pairA[3])		# 1-based in SAM
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
					if debug:
						agBAMDebug.debugger.debug("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d"
												  %(readsID,
												  contigA+":"+str(leftMostPosA),
												  contigB+":"+str(leftMostPosB),
												  int(alnLenA), int(alnLenB),
												  mapQA, mapQB, senseA, senseB, readLenA, readLenB))

					fracOvlA = alnLenA/readLenA
					fracOvlB = alnLenB/readLenB
					fracMismatchA = nMismatchesA/alnLenA
					fracMismatchB = nMismatchesB/alnLenB
					if (min(fracOvlA, fracOvlB) >= minFracOvl and				# minimum fraction of overlaps
						max(fracMismatchA, fracMismatchB) <= maxFracMismatch and	# maximum fraction of mismatches
						min(mapQA, mapQB) >= minMapQ):				# minimum mapping quality
						startA = leftMostPosA
						stopA = startA + int(alnLenA) - 1
						startB = leftMostPosB
						stopB = startB + int(alnLenB) - 1
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
					agBAMProgress.logger.info("%d parsed\t| %s\t| %.2f m" %(nReadsPairs, readsID, elapsedTime))
		with open(agBAMDone, 'w') as fDONE:
			fDONE.write("{}.bam done".format(agBAMOutbase))
	except KeyboardInterrupt:
		agBAMProgress.logger.info("Extract Joining-pairs INTERRUPTED by Keyboard")
		sys.exit(1)
	agBAMProgress.logger.info("{} reads pairs in {}".format(nReadsPairs, bamFile))
	agBAMProgress.logger.info("{} joining pairs parsed".format(nJoinPairs))
	agBAMProgress.logger.info("{} contig pairs given by these joining pairs".format(len(dContigPairs)))
	return (dContigPairs, nJoinPairs, bamFile)

def agouti_sam_main(bamFile, outDir, prefix,
					overwrite, minMapQ, minFracOvl,
					maxFracMismatch, nproc, debug=0):
	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_join_pairs")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)

	# 1. set up job and result queue
	config = moduleName, moduleOutDir, prefix, overwrite, minMapQ, minFracOvl, maxFracMismatch, debug
	print config
	jobs = []
	print bamFile
	for bam in bamFile:
		jobs.append( (bam, config) )

	def Start():
		print>>sys.stderr, 'Started a worker in %d from parent %d' %(os.getpid(), os.getppid())

	exec_pool = mp.Pool(nproc, initializer=Start)
	for i, res in enumerate(exec_pool.imap(worker, jobs)):
		#dContigPairs_each = collections.defaultdict(list)
		dContigPairs_each, nJoinPairs, bamFile = res
		print len(dContigPairs_each), nJoinPairs, bamFile
		if i == 0:
			dContigPairs = collections.defaultdict(list, dContigPairs_each)
		else:
			for k, v in dContigPairs_each.iteritems():
				dContigPairs[k].extend(v)
		print len(dContigPairs)
	#check_jp_entirity(dContigPairs, moduleOutDir)


	# 2. allocate BAMs in to job queue

	# 3. set up worker to read sam and support break and continue
	print "finished"
	sys.exit()

	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_join_pairs.progressMeter" %(prefix))
	agBAMOutAllJoinPairs = os.path.join(moduleOutDir, "%s.agouti.join_pairs.all.txt" %(prefix))
	agBAMProgress = agLOG.PROGRESS_METER(moduleName)
	if not os.path.exists(progressLogFile):
		agBAMProgress.add_file_handler(progressLogFile)
		agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")
	else:
		if not overwrite:
			agBAMProgress.add_file_handler(progressLogFile, 'a')
			dContigPairs = retrieve_joininng_pairs(agBAMProgress, agBAMOutAllJoinPairs)
			if dContigPairs is not None:
				return dContigPairs
			else:
				agBAMProgress.logger.info("Fail to pick up results from the previous run")
				agBAMProgress.logger.info("Re-processing the BAM file")
		else:
			agBAMProgress.add_file_handler(progressLogFile)
			agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")
			agBAMProgress.logger.info("Overwrite results from the previous run")

	agBAMDebug = None
	if debug:
		debugLogFile = os.path.join(moduleOutDir, "%s.agouti_join_pairs.debug" %(prefix))
		agBAMDebug = agLOG.DEBUG(moduleName, debugLogFile)

	# before running samtools, check its availability
	agBAMProgress.logger.info("check SAMtools")
	check_samtools(agBAMProgress)

	# runing samtools
	agBAMProgress.logger.info("run SAMtools")

	try:
		with open(agBAMOutAllJoinPairs, 'w') as fOUT:
			agBAMProgress.logger.info("# processed\t| Current Reads ID\t| Elapsed Time")
			if debug:
				agBAMDebug.debugger.debug("Reads_ID\tLocationA\tLocationB\tmapQA\tmapQB\tsenseA\tsenseB\treadLenA\treadLenB")
			startTime = time.time()
			dContigPairs = collections.defaultdict(list)
			nJoinPairs = 0
			nReadsPairs = 0
			for record in run_samtools(bamFile, agBAMProgress):
				tmpRecord = record.split("\n")
				pairA = tmpRecord[0].split("\t")
				pairB = tmpRecord[1].split("\t")
				readsID = pairA[0]
				contigA = pairA[2]
				contigB = pairB[2]
				mateCtgB = pairA[6]
				mateCtgA = pairB[6]
				nReadsPairs += 1
				# the first contidition makes sure
				# single end BAM are gonna have zero
				# joining-pairs extracted
				if contigA == "*" or contigB == "*":
					continue
				if pairA[0] == pairB[0] and contigA != contigB:
					alnLenA = getCIGAR(pairA[5])
					alnLenB = getCIGAR(pairB[5])
					leftMostPosA = int(pairA[3])		# 1-based in SAM
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
					if debug:
						agBAMDebug.debugger.debug("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d"
												  %(readsID,
												  contigA+":"+str(leftMostPosA),
												  contigB+":"+str(leftMostPosB),
												  int(alnLenA), int(alnLenB),
												  mapQA, mapQB, senseA, senseB, readLenA, readLenB))

					fracOvlA = alnLenA/readLenA
					fracOvlB = alnLenB/readLenB
					fracMismatchA = nMismatchesA/alnLenA
					fracMismatchB = nMismatchesB/alnLenB
					if (min(fracOvlA, fracOvlB) >= minFracOvl and				# minimum fraction of overlaps
						max(fracMismatchA, fracMismatchB) <= maxFracMismatch and	# maximum fraction of mismatches
						min(mapQA, mapQB) >= minMapQ):				# minimum mapping quality
						startA = leftMostPosA
						stopA = startA + int(alnLenA) - 1
						startB = leftMostPosB
						stopB = startB + int(alnLenB) - 1
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
					agBAMProgress.logger.info("%d parsed\t| %s\t| %.2f m" %(nReadsPairs, readsID, elapsedTime))
	except KeyboardInterrupt:
		agBAMProgress.logger.info("Extract Joining-pairs INTERRUPTED by Keyboard")
		sys.exit(1)

	agBAMProgress.logger.info("%d reads pairs in the give BAM" %(nReadsPairs))
	agBAMProgress.logger.info("%d joining pairs parsed" %(nJoinPairs))
	agBAMProgress.logger.info("%d contig pairs given by these joining pairs" %(len(dContigPairs)))
	if nJoinPairs == 0:
		agBAMProgress.logger.error("No joining pairs extracted")
		agBAMProgress.logger.error("Cannot SCAFFOLD without joining-pairs")
		sys.exit(1)
	else:
		agBAMProgress.logger.info("Succeeded")
	return dContigPairs

def get_joining_pairs(bamStream, outDir, prefix,
					  overwrite, minMapQ=5, minFracOvl=0.0,
					  maxFracMismatch=1.0, debug=0):

	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_join_pairs")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)

	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_join_pairs.progressMeter" %(prefix))
	agBAMOutAllJoinPairs = os.path.join(moduleOutDir, "%s.agouti.join_pairs.all.txt" %(prefix))
	agBAMProgress = agLOG.PROGRESS_METER(moduleName)
	if not os.path.exists(progressLogFile):
		agBAMProgress.add_file_handler(progressLogFile)
		agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")
	else:
		if not overwrite:
			agBAMProgress.add_file_handler(progressLogFile, 'a')
			dContigPairs = retrieve_joininng_pairs(agBAMProgress, agBAMOutAllJoinPairs)
			if dContigPairs is not None:
				return dContigPairs
			else:
				agBAMProgress.logger.info("Fail to pick up results from the previous run")
				agBAMProgress.logger.info("Re-processing the BAM file")
		else:
			agBAMProgress.add_file_handler(progressLogFile)
			agBAMProgress.logger.info("[BEGIN] Identifying joining pairs")
			agBAMProgress.logger.info("Overwrite results from the previous run")

	agBAMDebug = None
	if debug:
		debugLogFile = os.path.join(moduleOutDir, "%s.agouti_join_pairs.debug" %(prefix))
		agBAMDebug = agLOG.DEBUG(moduleName, debugLogFile)

	with open(agBAMOutAllJoinPairs, 'w') as fOUT:
		agBAMProgress.logger.info("# processed\t| Current Reads ID\t| Elapsed Time")
		if debug:
			agBAMDebug.debugger.debug("Reads_ID\tLocationA\tLocationB\tmapQA\tmapQB\tsenseA\tsenseB\treadLenA\treadLenB")
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
			readsID = pairA[0]
			contigA = pairA[2]
			contigB = pairB[2]
			nReadsPairs += 1
			if pairA[0] == pairB[0] and contigA != contigB:
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
				if debug:
					agBAMDebug.debugger.debug("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d" %(readsID,
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
				agBAMProgress.logger.info("%d parsed\t| %s\t| %.2f m" %(nReadsPairs, readsID, elapsedTime))

	agBAMProgress.logger.info("%d joining pairs parsed" %(nJoinPairs))
	agBAMProgress.logger.info("%d contig pairs given by these joining pairs" %(len(dContigPairs)))
	if nJoinPairs == 0:
		agBAMProgress.logger.error("No joining pairs extracted")
		agBAMProgress.logger.error("Cannot SCAFFOLD without joining-pairs")
		sys.exit(1)
	else:
		agBAMProgress.logger.info("Succeeded")
	return dContigPairs
