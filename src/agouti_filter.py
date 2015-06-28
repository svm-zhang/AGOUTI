import sys
import os
import collections

from lib import agouti_gff as agGFF

def createNewGenes(geneModels, index, ctg, start, stop):		# index is 0-based
	geneModel = agGFF.AGOUTI_GFF()
	geneModel.setGene("AGOUTI.%s.%d" %(ctg, len(geneModels)),
					   -1, -1, 0)
	geneModel.lcds = [start, stop]
	geneModel.setContigID = ctg
	geneModel.setProgram("AGOUTI")
	if index == 0:
		geneModels =  [geneModel] + geneModels
	elif index == len(geneModels)-1:
		geneModels = geneModels.append(geneModel)
	else:
		geneModels = geneModels[0:index] + [geneModel] + geneModels[index:]
#	print geneModels[index].geneID, geneModels[index].lcds, ctg
	return geneModels

def matchGene(geneModels, start, stop):
	ctg = geneModels[0].ctgID
	if len(geneModels) == 1:
		# only one gene on this contig
		# interval hits the gene
		if start >= geneModels[0].geneStart and stop <= geneModels[0].geneStop:
			print "only gene hit", ctg, geneModels[0].geneStart, geneModels[0].geneStop, start, stop
			return geneModels, 0
#			return geneModels[0], 0
		# interval spans the gene
		elif ((start < geneModels[0].geneStart and stop > geneModels[0].geneStart) or
			  (start >= geneModels[0].geneStart and start <= geneModels[0].geneStop and stop > geneModels[0].geneStop) or
			  (start < geneModels[0].geneStart and stop >= geneModels[0].geneStart and stop <= geneModels[0].geneStop)):
			print "only gene span", ctg, geneModels[0].geneStart, geneModels[0].geneStop, start, stop
			return geneModels, 0
#			return geneModels[0], 0
		# outside the gene, require create new gene model
		else:
			if stop < geneModels[0].geneStart:
				createNewGenes(geneModels, 0, ctg, start, stop)
				print "only gene before", ctg, geneModels[0].geneStart, geneModels[0].geneStop, start, stop
				return geneModels, 0
#				return geneModels[0], 0
			elif start > geneModels[0].geneStop:
				createNewGenes(geneModels, len(geneModels)-1, ctg, start, stop)
				print "only gene after", ctg, geneModels[0].geneStart, geneModels[0].geneStop, start, stop
				return geneModels, len(geneModels)-1
#				return geneModels[0], len(geneModels[0])-1
	for i in xrange(0, len(geneModels)):
		curGeneModel = geneModels[i]
		# first gene
		if i == 0:
			# interval hits the gene
			if start >= curGeneModel.geneStart and stop <= curGeneModel.geneStop:
				print "first hit", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
				return geneModels, i
#				return curGeneModel, i
			# interval spans the gene
			elif ((start < curGeneModel.geneStart and stop > curGeneModel.geneStart) or
				  (start >= curGeneModel.geneStart and start <= curGeneModel.geneStop and stop > curGeneModel.geneStop) or
				  (start < curGeneModel.geneStart and stop >= curGeneModel.geneStart and stop <= curGeneModel.geneStop)):
				print "first span", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
				return geneModels, i
#				return curGeneModel, i
			else:
				# interval before the first gene
				if stop < curGeneModel.geneStart:
					print "first before", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
#					geneModels = createNewGenes(geneModels, 0, ctg, start, stop)
					return geneModels, 0
#					return curGeneModel, 0
		# last gene
		elif i == len(geneModels)-1:
			# interval hits the gene
			if start >= curGeneModel.geneStart and stop <= curGeneModel.geneStop:
				print "last hit", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
				return geneModels, i
#				return curGeneModel, i
			# interval spans the gene
			elif ((start < curGeneModel.geneStart and stop > curGeneModel.geneStart) or
				  (start >= curGeneModel.geneStart and start <= curGeneModel.geneStop and stop > curGeneModel.geneStop) or
				  (start < curGeneModel.geneStart and stop >= curGeneModel.geneStart and stop <= curGeneModel.geneStop)):
				print "last span", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
				return geneModels, i
#				return curGeneModel, i
			else:
				# interval after the last gene
				if start > curGeneModel.geneStop:
					print "last after", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
#					geneModels = createNewGenes(geneModels, len(geneModels)-1, ctg, start, stop)
					return geneModels, i
#					return curGeneModel, -1
				# intergenic region before the last gene
				elif stop <= curGeneModel.geneStart and start >= preGeneModel.geneStop:
					print "last before", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
#					geneModels = createNewGenes(geneModels, len(geneModels)-2, ctg, start, stop)
					return geneModels, i-1
#					return curGeneModel, i-1
		else:
			# interval hits the gene
			if start >= curGeneModel.geneStart and stop <= curGeneModel.geneStop:
				print "internal hit", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
				return geneModels, i
#				return curGeneModel, i
			# interval spans the gene
			elif ((start < curGeneModel.geneStart and stop > curGeneModel.geneStart) or
				  (start >= curGeneModel.geneStart and start <= curGeneModel.geneStop and stop > curGeneModel.geneStop) or
				  (start < curGeneModel.geneStart and stop >= curGeneModel.geneStart and stop <= curGeneModel.geneStop)):
				print "internal span", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
				return geneModels, i
#				return curGeneModel, i
			# interval in between two genes
			else:
				if stop < curGeneModel.geneStart and start > preGeneModel.geneStop:
					print "intergenic", ctg, geneModels[i].geneStart, geneModels[i].geneStop, start, stop
#					geneModels = createNewGenes(geneModels, i, ctg, start, stop)
					return geneModels, i
#					return curGeneModel, i
		preGeneModel = curGeneModel
#	all of the above conditions should include every cases,
#	alway return this if there is a bug
	return geneModels, -2			# -2 encodes cases that failed to consider
#	return None, -2

def recordGeneIndex(contigs, geneModels):
	for i in range(len(contigs)):
		contig = contigs[i]
		geneModel = geneModels[i]
		geneID = geneModel.geneID

def cleanContigPairs(dContigPairs, dGFFs, joinPairsFile):
	sys.stderr.write("Filtering joining pairs ... \n")
	dCtgPair2GenePair = collections.defaultdict()
	dMappedPos = collections.defaultdict()
	daddedModels = collections.defaultdict(list)
	fOUT = open(joinPairsFile, 'w')
	for ctgPair, pairInfo in dContigPairs.items():
		ctgA = ctgPair[0]
		ctgB = ctgPair[1]
		print ">%s %s" %(ctgA, ctgB)
		pairToRemove = []
		for i in xrange(len(pairInfo)):
			startA = pairInfo[i][0]
			stopA = pairInfo[i][2]
			startB = pairInfo[i][1]
			stopB = pairInfo[i][3]
			senseA = pairInfo[i][4]
			senseB = pairInfo[i][5]
			readID = pairInfo[i][-1]
			geneIndexA, geneIndexB = -1, -1
			if ctgA in dGFFs:
				dGFFs[ctgA], geneIndexA = matchGene(dGFFs[ctgA], startA, stopA)
#				geneModelA, geneIndexA = matchGene(dGFFs[ctgA], startA, stopA)
			else:
				dGFFs[ctgA] = createNewGenes(dGFFs[ctgA], 0, ctgA, startA, stopA)
				geneIndexA = 0
			if ctgB in dGFFs:
				dGFFs[ctgB], geneIndexB = matchGene(dGFFs[ctgB], startB, stopB)
#				geneModelB, geneIndexB = matchGene(dGFFs[ctgB], startB, stopB)
			else:
				dGFFs[ctgB] = createNewGenes(dGFFs[ctgB], 0, ctgB, startB, stopB)
				geneIndexB = 0

			print geneIndexA, geneIndexB
			geneModelA = dGFFs[ctgA][geneIndexA]
			geneModelB = dGFFs[ctgB][geneIndexB]
			nGeneToLeftA, nGeneToRightA = -1, -1
			nGeneToLeftB, nGeneToRightB = -1, -1
			# use geneIndex as check conditions
			if geneIndexA != -1 and geneIndexB != -1:
				nGeneToLeftA = geneIndexA
				nGeneToRightA = len(dGFFs[ctgA]) - 1 - geneIndexA
				nGeneToLeftB = geneIndexB
				nGeneToRightB = len(dGFFs[ctgB]) - 1 - geneIndexB
				print pairInfo[0], geneIndexA, nGeneToLeftA, nGeneToRightA, len(dGFFs[ctgA]), geneIndexB, nGeneToLeftB, nGeneToRightB, len(dGFFs[ctgB])
				# case where RR with 55
				if (senseA == senseB and senseA == '-' and
					nGeneToLeftA <= nGeneToRightA and
					nGeneToLeftB <= nGeneToRightB and
					nGeneToLeftA + nGeneToLeftB <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, geneIndexA, len(dGFFs[ctgA]), geneIndexB, len(dGFFs[ctgB]), nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
					if (ctgA, ctgB) not in dCtgPair2GenePair:
						dCtgPair2GenePair[ctgA, ctgB] = [geneModelA, geneModelB]
#					dMappedPos = recordGeneIndex(list(ctgA, ctgB), list(geneModelA, geneModelB))
				elif (senseA == senseB and senseA == '+' and
					nGeneToLeftA >= nGeneToRightA and
					nGeneToLeftB >= nGeneToRightB and
					nGeneToRightA + nGeneToRightB <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, geneIndexA, len(dGFFs[ctgA]), geneIndexB, len(dGFFs[ctgB]), nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
					if (ctgA, ctgB) not in dCtgPair2GenePair:
						dCtgPair2GenePair[ctgA, ctgB] = [geneModelA, geneModelB]
				elif (senseA == '+' and senseB == '-' and
					nGeneToLeftA >= nGeneToRightA and
					nGeneToLeftB <= nGeneToRightB and
					nGeneToRightA + nGeneToLeftB <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
					if (ctgA, ctgB) not in dCtgPair2GenePair:
						dCtgPair2GenePair[ctgA, ctgB] = [geneModelA, geneModelB]
				elif (senseA == '-' and senseB == '+' and
					nGeneToLeftA <= nGeneToRightA and
					nGeneToLeftB >= nGeneToRightB and
					nGeneToRightB + nGeneToLeftA <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
					if (ctgA, ctgB) not in dCtgPair2GenePair:
						dCtgPair2GenePair[ctgA, ctgB] = [geneModelA, geneModelB]
				else:
					print "\t", "both else", pairInfo[i], geneIndexA, geneIndexB, nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
					print "\t", "both else", nGeneToLeftA + nGeneToRightB, nGeneToLeftA <= nGeneToRightA, nGeneToLeftB >= nGeneToRightB
					pairToRemove.append(i)
			elif geneIndexA == -1 and geneIndexB == -1:
				# this case should not happen
				print "\t", "None", pairInfo[i]
				print "\t", "None", ctgA, ctgB, "None", "None", "N/A", "N/A", "N/A", "N/A", senseA, senseB
				fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
			elif geneIndexA == -2 or geneIndexB == -2:
				sys.stderr.write("%s\n" %(str(pairInfo[i])))
				sys.stderr.write("ding ding ding! check is in emergency!")
				sys.exit(1)
		if len(pairInfo) - len(pairToRemove) < 3:			# 5 (number of links) is a magical number, provide an argument in command line
			del dContigPairs[ctgPair]
			print "\t", "remove enitre pair"
		else:
			print "### removed links: ", len(pairToRemove)
			print "### links in total: ", len(pairInfo), "\n"
			tmp = []
			for i in xrange(len(pairInfo)):
				if i not in pairToRemove:
					tmp.append(pairInfo[i])
				else:
					pairToRemove.remove(i)
			dContigPairs[ctgPair] = tmp
	fOUT.close()
	return dCtgPair2GenePair
