import sys
import os
import collections
import operator
import logging

from lib import agouti_gff as agGFF
from lib import agouti_log as agLOG

def merge_intervals(intervals):
	if len(intervals) == 0:
		return []
	mergedIntervals = []
	sortIntervals = sorted(intervals)
	for higher in sortIntervals:
		if not mergedIntervals:
			mergedIntervals.append(higher)
		else:
			lower = mergedIntervals[-1]
			if higher[0] <= lower[1]:
				upper_bound = max(lower[1], higher[1])
				mergedIntervals[-1] = (lower[0], upper_bound)
			else:
				mergedIntervals.append(higher)
	print "\t    merged intervals", mergedIntervals
	return mergedIntervals

def create_fake_genes(geneModels, index, ctg, intervals):
	geneModel = agGFF.AGOUTI_GFF()
	fake = 1
	start = intervals[0][0]
	stop = intervals[-1][1]
	geneModel.setContigID = ctg
	geneModel.setGene("%s_AGOUTI_fake_%d" %(ctg, len(geneModels)),
					  start, stop, fake)
	geneModel.setProgram("AGOUTI")
	geneModel.lcds = [pos for interval in intervals for pos in interval]
	if index == 0 or index == -1:
		geneModels =  [geneModel] + geneModels
	elif index == len(geneModels):
		geneModels = geneModels + [geneModel]
	else:
		geneModels = geneModels[0:index] + [geneModel] + geneModels[index:]
	print "\tcreate fake gene", geneModel.geneID, geneModel.geneStart, geneModel.geneStop, geneModel.lcds
	return geneModels

def find_interval_belongings(intervalQry, mergedIntervals):
	for interval in mergedIntervals:
		if intervalQry[1] >= interval[0] and interval[1] >= intervalQry[0]:
			return interval

def find_overlap(interval1, interval2):
	if interval1[1] >= interval2[0] and interval2[1] >= interval1[0]:
		# 0: two intervals overlap
		return 0
	elif interval1[1] < interval2[0]:
		# -1: interval1 on the left side of interval2
		return -1
	elif interval1[0] > interval2[1]:
		# 1: interval1 on the right side of interval2
		return 1

def find_gene_overlap(interval, geneModel_5, geneModel_3):
	print "find gene overlap for interval", interval
	if geneModel_5.geneID == geneModel_3.geneID:
		print "single gene contig", geneModel_5.geneID, geneModel_5.geneStart, geneModel_5.geneStop
		ovl = find_overlap(interval, (geneModel_5.geneStart, geneModel_5.geneStop))
		print "ovl", ovl
		if ovl == 0:
			# relatviePos, geneIndex, end
			return 0, 0
		elif ovl == -1:
			return -1, 0
		elif ovl == 1:
			return 1, 0
	else:
		print "multiple gene contig"
		print "\tgeneModel_5", geneModel_5.geneID, geneModel_5.geneStart, geneModel_5.geneStop
		print "\tgeneModel_3", geneModel_3.geneID, geneModel_3.geneStart, geneModel_3.geneStop
		ovl_5 = find_overlap(interval, (geneModel_5.geneStart, geneModel_5.geneStop))
		ovl_3 = find_overlap(interval, (geneModel_3.geneStart, geneModel_3.geneStop))
		print "ovl_5", ovl_5
		print "ovl_3", ovl_3
		if ovl_5 == 0 and ovl_3 != 0:
			return 0, 5
		elif ovl_5 != 0 and ovl_3 == 0:
			return 0, 3
		else:
			if ovl_5 == -1 and ovl_3 == -1:
				return -1, 5
			elif ovl_5 == 1 and ovl_3 == 1:
				return 1, 3
			elif ovl_5 == 1 and ovl_3 == -1:
				return -2, -2
			elif ovl_5 == 0 and ovl_3 == 0:
				#!!! should merge two genes
				return 0, 0

def mapping_to_geneModel(geneModelsA, geneModelsB,
						 mapIntervalsA, mapIntervalsB,
						 mergedIntervalsA, mergedIntervalsB,
						 senses):
	geneModelA_5 = geneModelsA[0]
	geneModelA_3 = geneModelsA[-1]
	geneModelB_5 = geneModelsB[0]
	geneModelB_3 = geneModelsB[-1]
	dMappingSupport = {}
	dSenses = {}
	print "mapping_to_geneModel all senses", collections.Counter(senses)
	for i in xrange(len(mapIntervalsA)):
		intervalA = mapIntervalsA[i]
		mergedIntervalA = find_interval_belongings(intervalA, mergedIntervalsA)
		intervalB = mapIntervalsB[i]
		mergedIntervalB = find_interval_belongings(intervalB, mergedIntervalsB)
		if (mergedIntervalA, mergedIntervalB) not in dMappingSupport:
			dMappingSupport[mergedIntervalA, mergedIntervalB] = 1
			dSenses[mergedIntervalA, mergedIntervalB] = [senses[i]]
		else:
			dMappingSupport[mergedIntervalA, mergedIntervalB] += 1
			dSenses[mergedIntervalA, mergedIntervalB] += [senses[i]]

	createOnBoth = []
	dIndex2Cluster = collections.defaultdict(list)
	for k in sorted(dMappingSupport.items(), key=operator.itemgetter(1), reverse=True):
		clusterA = k[0][0]
		clusterB = k[0][1]
		nSupport = k[1]
		print clusterA, clusterB, nSupport
		print "A"
		geneIndexA, endA = find_gene_overlap(clusterA, geneModelA_5, geneModelA_3)
		print "B"
		geneIndexB, endB = find_gene_overlap(clusterB, geneModelB_5, geneModelB_3)
		print "geneIndexA", geneIndexA, endA, "geneIndexB", geneIndexB, endB
		if geneIndexA == 0 and geneIndexB == 0:
			# save and go on
			return (geneIndexA, geneIndexB, endA, endB, [], [], dSenses[k[0]])
		elif geneIndexA == -2 or geneIndexB == -2:
			return None
		elif geneIndexA != 0 and geneIndexB == 0:
			return (geneIndexA, geneIndexB, endA, endB, [clusterA], [], dSenses[k[0]])
		elif geneIndexA == 0 and geneIndexB != 0:
			return (geneIndexA, geneIndexB, endA, endB, [], [clusterB], dSenses[k[0]])
		elif geneIndexA != 0 and geneIndexB != 0:
			createOnBoth.append((clusterA, clusterB))
			dIndex2Cluster[geneIndexA, geneIndexB, endA, endB].append((clusterA, clusterB))

	print "createOnBoth", createOnBoth
	print "dIndex2Cluster", dIndex2Cluster
	if len(createOnBoth) > 0:
		# create gene on both contig  using all mapped intervals
		print "createonboth"
		if len(dIndex2Cluster) == 1:
			genePair = list(dIndex2Cluster.keys()[0])
			intervalsA = merge_intervals([cluster[0] for x in dIndex2Cluster.itervalues() for cluster in x])
			intervalsB = merge_intervals([cluster[1] for x in dIndex2Cluster.itervalues() for cluster in x])
			tmp = []
			for k, v in dIndex2Cluster.iteritems():
				print "k", k, "v", v
				for i in range(len(v)):
					tmp += dSenses[v[i]]
			print "intervalsA", intervalsA
			print "intervalsB", intervalsB
			print "tmp", tmp
			return tuple(genePair + [intervalsA] + [intervalsB] + [tmp])
		else:
			print "index conflict"

def get_genePair_for_contigPair(dGFFs, ctgA, ctgB, mapIntervalsA, mapIntervalsB, senses):
	geneModelsA = []
	geneModelsB = []
	createA, createB = 0, 0
	mergedIntervalsA = merge_intervals(mapIntervalsA)
	if ctgA not in dGFFs:
		# create genes on contig A using all intervals
		print "create gene on A"
		dGFFs[ctgA] = create_fake_genes([], 0, ctgA, mergedIntervalsA)
		createA = 1
	mergedIntervalsB = merge_intervals(mapIntervalsB)
	if ctgB not in dGFFs:
		# create genes on contig B using all intervals
		print "create gene on B"
		dGFFs[ctgB] = create_fake_genes([], 0, ctgB, mergedIntervalsB)
		createB = 1
	print "createA", createA, "createB", createB
	if createA and createB:
		# no need to map to gene model, return gene pair directly
		print "create both"
		return (0, 0, 0, 0, [], [], senses)
	geneModelsA = dGFFs[ctgA]
	geneModelsB = dGFFs[ctgB]
	print "get_genePair_for_contigPair all senses", senses
	genePair = mapping_to_geneModel(geneModelsA, geneModelsB,
									mapIntervalsA, mapIntervalsB,
									mergedIntervalsA, mergedIntervalsB,
									senses)
	print "contig Pair", ctgA, ctgB, "genePair", genePair
	return genePair

def set_module_name(name):
	global moduleName
	moduleName = name

def enforce_filters(dContigPairs, dGFFs,
					moduleOutDir, prefix, minSupport):
	moduleProgressLogFile = os.path.join(moduleOutDir, "%s.agouti_filter.progressMeter" %(prefix))
	moduleDebugLogFile = os.path.join(moduleOutDir, "%s.agouti_filter.debug" %(prefix))
	moduleOutputFile = os.path.join(moduleOutDir, "%s.agouti.join_pairs.filtered.txt" %(prefix))
	moduleProgressLogger = agLOG.AGOUTI_LOG(moduleName).create_logger(moduleProgressLogFile)
	moduleDEBUGLogger = agLOG.AGOUTI_DEBUG_LOG(moduleName+"_DEBUG").create_logger(moduleDebugLogFile)

	moduleProgressLogger.info("[BEGIN] Filtering joining pairs")
	dCtgPair2GenePair = collections.defaultdict()
	dMappedPos = collections.defaultdict()
	daddedModels = collections.defaultdict(list)
	nFail4Combination = 0
	nFailGeneModel = 0
	fOUT = open(moduleOutputFile, 'w')
	for ctgPair, pairInfo in dContigPairs.items():
		if len(pairInfo) < minSupport:
			del dContigPairs[ctgPair]
			continue
		ctgA = ctgPair[0]
		ctgB = ctgPair[1]
		print ">%s %s" %(ctgA, ctgB)
		pairToRemove = []
		mapIntervalsA = []
		mapIntervalsB = []
		pairs = []
		senses = []
		keep = 0
		for i in xrange(len(pairInfo)):
			startA, startB, stopA, stopB, senseA, senseB, readID = pairInfo[i]
			mapIntervalsA += [(startA, stopA)]
			mapIntervalsB += [(startB, stopB)]
			pairs += [(startA, stopA, startB, stopB)]
			senses += [(senseA, senseB)]
		genePair = get_genePair_for_contigPair(dGFFs, ctgA, ctgB, mapIntervalsA, mapIntervalsB, senses)
		geneModelsA = dGFFs[ctgA]
		geneModelsB = dGFFs[ctgB]
		if genePair is not None:
			geneIndexA, geneIndexB, endA, endB, intervalsA, intervalsB, senses = genePair
			sensesCounter = collections.Counter(senses)
			print "sensesCounter", sensesCounter
			if geneIndexB != 0:
				# create gene model according to endB using intervalsB
				if geneIndexB == -1 and (endB == 5 or endB == 0):
					dGFFs[ctgB] = create_fake_genes(geneModelsB, 0, ctgB, intervalsB)
					geneIndexB = 0
					endB = 5
				elif geneIndexB == 1 and (endB == 3 or endB == 0):
					dGFFs[ctgB] = create_fake_genes(geneModelsB, len(geneModelsB), ctgB, intervalsB)
					geneIndexB = len(dGFFs[ctgB]) - 1
					endB = 3
			else:
				if endB == 0:
					endB = 5
				elif endB == 3:
					geneIndexB = len(dGFFs[ctgB])-1
			if geneIndexA != 0:
				# create gene model according to endA using intervalsA
				if geneIndexA == -1 and (endA == 5 or endA == 0):
					dGFFs[ctgA] = create_fake_genes(geneModelsA, 0, ctgA, intervalsA)
					geneIndexA = 0
					endA = 5
				elif geneIndexA == 1 and (endA == 3 or endA == 0):
					dGFFs[ctgA] = create_fake_genes(geneModelsA, len(geneModelsA), ctgA, intervalsA)
					geneIndexA = len(dGFFs[ctgA]) - 1
					endA = 3
			else:
				if endA == 0:
					endA = 3
				elif endA == 3:
					geneIndexA = len(dGFFs[ctgA])-1
			print "contig Pair", ctgA, ctgB, "genePair", genePair
			print "# model on ctgA", len(dGFFs[ctgA]), "# model on ctgB", len(dGFFs[ctgB])
			print "geneIndexA", geneIndexA, "endA", endA, "geneIndexB", geneIndexB, "endB", endB
			sense = sorted(sensesCounter.items(), key=operator.itemgetter(1), reverse=True)[0][0]
			print "sense", sense
			if (geneIndexA == len(dGFFs[ctgA])-1 and endA == 3) and \
			   (geneIndexB == 0 and endB == 5) and sense == ('+', '-'):
					# FR + 3'-5'
					keep = 1
			elif (geneIndexA == 0 and endA == 5) and \
				 (geneIndexB == 0 and endB == 5) and sense == ('-', '-'):
					# RR + 5'-5'
					keep = 1
			elif (geneIndexA == len(dGFFs[ctgA])-1 and endA == 3) and \
				 (geneIndexB == len(dGFFs[ctgB])-1 and endB == 3) and \
				 sense == ('+', '+'):
					# FF + 3'-3'
					keep = 1
			elif (geneIndexA == 0 and endA == 5) and \
				 (geneIndexB == len(dGFFs[ctgB])-1 and endB == 3) and \
				 sense == ('-', '+'):
					# RF + 5'-3'
					keep = 1
			elif (geneIndexA == 0 and (endA == 0 or endA == 3)) and \
				 (geneIndexB == 0 and (endB == 0 or endB == 5)) and \
				 sense == ('+', '-'):
					# only one gene on the contig
					# it doesn't matter which end
					keep = 1
			print "****keep", keep
			if keep:
				geneA = dGFFs[ctgA][geneIndexA]
				geneB = dGFFs[ctgB][geneIndexB]
				dCtgPair2GenePair[ctgA, ctgB] = [geneA, geneB]
				print "geneA", geneA.geneID, geneA.geneStart, geneA.geneStop
				print "geneB", geneB.geneID, geneB.geneStart, geneB.geneStop
				print "sense", sense
				senseA = sense[0]
				senseB = sense[1]
				for i in xrange(len(pairInfo)):
					startA, startB, stopA, stopB, _, _, readID = pairInfo[i]
					intervalA = (startA, stopA)
					intervalB = (startB, stopB)
					print "intervalA", intervalA, "intervalB", intervalB
					if len(intervalsA) == 0:
						if len(intervalsB) == 0:
							print "use all"
							fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
						else:
							print "use all A, not all B"
							overlap = find_overlap(intervalB, (geneB.geneStart, geneB.geneStop))
							if overlap == 0:
								fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
					else:
						if len(intervalsB) == 0:
							print "use all B, not all A"
							overlap = find_overlap(intervalA, (geneA.geneStart, geneA.geneStop))
							if overlap == 0:
								fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
						else:
							print "not all Both"
							overlapA = find_overlap(intervalA, (geneA.geneStart, geneA.geneStop))
							overlapB = find_overlap(intervalB, (geneB.geneStart, geneB.geneStop))
							if overlapA == 0 and overlapB == 0:
								fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
			else:
				nFail4Combination += 1
#			if len(sensesCounter) == 1:
#				sense = sensesCounter.keys()[0]
#			else:
#				print "multiple sense pairs"
#				senses = sorted(sensesCounter.items(), key=operator.itemgetter(1), reverse=True)[0:2]
#				print "senses", senses
#				ratio = float(senses[0][1])/(senses[0][1]+senses[1][1])
#				print "ratio", ratio
		else:
			nFailGeneModel += 1
			print "genePair is None"
	fOUT.close()
	moduleProgressLogger.info("%d contig pairs filtered for spanning across >1 gene models" %(nFailGeneModel))
	moduleProgressLogger.info("%d contig pairs filtered for not being one of the four combinations" %(nFail4Combination))
	moduleProgressLogger.info("Succeeded")
	return dCtgPair2GenePair, moduleOutputFile

