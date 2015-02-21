#!/usr/bin/python

import sys
import os
import optparse
import string
import collections
import math
from numpy import zeros, int16

parpath = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), os.pardir))
if os.path.exists(os.path.join(parpath, "lib", "__init__.py")):
	sys.path.insert(0, parpath)

from lib import agouti_gff as agff
from lib import agouti_sam as asam

versionString = "RNAPATH*" # modified from RNAPATH by luting zhuo

def compNT(nt):
	""" returns the complementary basepair to base nt
	"""
	compDict = { "A": "T",
			 "T": "A",
			 "G": "C",
			 "C": "G",
			 "S": "S",
			 "W": "W",
			 "R": "Y",
			 "Y": "R",
			 "M": "K",
			 "K": "M",
			 "H": "D",
			 "D": "H",
			 "B": "V",
			 "V": "B",
			 "N": "N",
			 "a": "t",
			 "t": "a",
			 "g": "c",
			 "c": "g",
			 "n": "n",
			 "z": "z"
	}

	return compDict.get(nt, "N")

def complement(sequence, length=-1):
	""" returns the complement of the sequence.
	"""
	newSeq = ""
	seqLength = len(sequence)
	if length == seqLength or length < 0:
		seqList = list(sequence)
		seqList.reverse()
		return "".join(map(compNT, seqList))

#TODO: this seems to want to deal with case where length is more than
# sequence length except that a negative index on a sequence is fine
# index will only be overrun if length is negative but that case is
# handled above
	for index in range(seqLength - 1,seqLength - length - 1, -1):
		try:
			newSeq += compNT(sequence[index])
		except:
			newSeq += "N"

	return newSeq

def main(argv=None):
	if not argv:
		argv = sys.argv

	usage = "python %prog incontigfile distalPairs outpathfile outcontigfile [--prefix string] [--overlap bp]"

	if len(argv) < 5:
		print usage
		sys.exit(0)

	icontig = argv[1]
#	distalPairsfile = argv[2]
	isam = argv[2]
	igff = argv[3]
	outpathfilename = argv[4]
	outcontigfilename = argv[5]

	rnaPath(icontig, isam, igff, outpathfilename,
			outcontigfilename)

def gff2GeneModels(igff):
	dGFFs = collections.defaultdict(list)
	nGene = 0
	with open(igff, 'r') as fIN:
		for line in fIN:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				if tmp_line[2] == "gene":
					nGene += 1
		lobj_GeneModels = [agff.AGOUTI_GFF() for i in xrange(nGene)]
		geneIndex = -1
		fIN.seek(0)
		for line in fIN:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				if tmp_line[2] == "gene":
					geneIndex += 1
					if geneIndex == 0:
						lobj_GeneModels[geneIndex].setGene(tmp_line[8].split('=')[1],
														   int(tmp_line[3]),
														   int(tmp_line[4]), 0)
					else:
						preCtgID = lobj_GeneModels[geneIndex-1].ctgID
						preGeneID = lobj_GeneModels[geneIndex-1].geneID
						dGFFs[preCtgID].append(lobj_GeneModels[geneIndex-1])
						lobj_GeneModels[geneIndex].setGene(tmp_line[8].split('=')[1],
														   int(tmp_line[3]),
														   int(tmp_line[4]), 0)
					lobj_GeneModels[geneIndex].setProgram(tmp_line[1])
					lobj_GeneModels[geneIndex].setContigID(tmp_line[0])
					lobj_GeneModels[geneIndex].setStrand(tmp_line[6])
				elif tmp_line[2] == "stop_codon":
					lobj_GeneModels[geneIndex].setStopCodon()
				elif tmp_line[2] == "start_codon":
					lobj_GeneModels[geneIndex].setStartCodon()
				elif tmp_line[2] == "CDS":
					lobj_GeneModels[geneIndex].updateCDS(int(tmp_line[3]), int(tmp_line[4]))
		dGFFs[lobj_GeneModels[geneIndex].ctgID].append(lobj_GeneModels[geneIndex])

	nGeneModels = 0
	for k, v in sorted(dGFFs.items()):
		nGeneModels += len(v)
	sys.stderr.write("Number of Gene Models parsed: %d\n" %(nGeneModels))
	return dGFFs

def addGeneModel(ctg, daddedModels, start, stop):
	if ctg not in daddedModels:
		geneModel = agff.AGOUTI_GFF()
		geneModel.lcds = [start, stop]
		geneModel.setContigID = ctg
		geneModel.setProgram("AGOUTI")
		daddedModels[ctg] = [geneModel]
		return geneModel, daddedModels
	else:
		for i in range(len(daddedModels[ctg])):
			tmpModel = daddedModels[ctg][i]
			print tmpModel.lcds
			# make another gene
			if ((start < stop and tmpModel.lcds[0] - stop >= 10000) or
				(start < stop and start - tmpModel.lcds[-1] >= 10000)):
				geneModel = agff.AGOUTI_GFF()
				geneModel.lcds = [start, stop]
				geneModel.setContigID = ctg
				geneModel.setProgram("AGOUTI")
				daddedModels[ctg].append(geneModel)
				return geneModel, daddedModels
			for j in range(0, len(tmpModel.lcds), 2):
				print tmpModel.lcds[j], tmpModel.lcds[j+1]
				# outside exon
				if (start < stop and tmpModel.lcds[j] > stop):
					# add before the first exon
					if j == 0:
						tmpModel.lcds = [start, stop] + tmpModel.lcds
						return tmpModel, daddedModels
					# squeeze in a new exon
					else:
						tmpModel.lcds = tmpModel.lcds[:j] + [start, stop] + tmpModel.lcds[j:]
						return tmpModel, daddedModels
				# 3' to the last exon
				elif (j == len(tmpModel.lcds)-2 and start < stop and start > tmpModel.lcds[j+1]):
					tmpModel.lcds += [start, stop]
					return tmpModel, daddedModels
				# overlap exon
				elif (start < tmpModel.lcds[j] and stop <= tmpModel.lcds[j+1] and stop > tmpModel.lcds[j]):
					tmpModel.lcds[j] = start
					return tmpModel, daddedModels
				# overlap exon
				elif (start < tmpModel.lcds[j+1] and start >= tmpModel.lcds[j] and stop > tmpModel.lcds[j+1]):
					tmpModel.lcds[j+1] = stop
					return tmpModel, daddedModels
				# within exon
				elif (start >= tmpModel.lcds[j] and stop <= tmpModel.lcds[j+1]):
					return tmpModel, daddedModels
				# span exon
				elif (start < tmpModel.lcds[j] and stop > tmpModel.lcds[j+1]):
					tmpModel.lcds[j] = start
					tmpModel.lcds[j+1] = stop
					return tmpModel, daddedModels
				else:
					if j == len(tmpModel.lcds)-2:
						print "This should not happen; debug!", ctg, start, stop

def createNewGenes(geneModels, index, ctg, start, stop):		# index is 0-based
	geneModel = agff.AGOUTI_GFF()
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

def matchGene(geneModels, start, stop):
	ctg = geneModels[0].ctgID
	# only one gene on this contig
	# interval hits the gene
	if start >= geneModels[0].geneStart and stop <= geneModels[0].geneStop:
		return 0
#		return geneModels[0], 0
	# interval spans the gene
	elif ((start < geneModels[0].geneStart and stop > geneModels[0].geneStart) or
		  (start >= geneModels[0].geneStart and stop > geneModels[0].geneStop) or
		  (start < geneModels[0].geneStart and stop <= geneModels[0].geneStop)):
		return 0
#		return geneModels[0], 0
	# outside the gene, require create new gene model
	else:
		if stop < geneModels[0].geneStart:
			createNewGenes(geneModels, 0, ctg, start, stop)
			return 0
#			return geneModels[0], 0
		elif start > geneModels[0].geneStop:
			createNewGenes(geneModels, len(geneModels)-1, ctg, start, stop)
			return len(geneModels)-1
#			return geneModels[0], len(geneModels[0])-1
	for i in xrange(0, len(geneModels)):
		curGeneModel = geneModels[i]
#		if i == 0 and i == len(geneModels) - 1:
		# first gene
		if i == 0:
			# interval hits the gene
			if start >= curGeneModel.geneStart and stop <= curGeneModel.geneStop:
				return i
#				return curGeneModel, i
			# interval spans the gene
			elif ((start < curGeneModel.geneStart and stop > curGeneModel.geneStart) or
				  (start >= curGeneModel.geneStart and stop > curGeneModel.geneStop) or
				  (start < curGeneModel.geneStart and stop <= curGeneModel.geneStop)):
				return i
#				return curGeneModel, i
			else:
				# interval before the first gene
				if stop < curGeneModel.geneStart:
					createNewGenes(geneModels, 0, ctg, start, stop)
					return 0
#					return curGeneModel, 0
		# last gene
		elif i == len(geneModels)-1:
			# interval hits the gene
			if start >= curGeneModel.geneStart and stop <= curGeneModel.geneStop:
				return i
#				return curGeneModel, i
			# interval spans the gene
			elif ((start < curGeneModel.geneStart and stop > curGeneModel.geneStart) or
				  (start >= curGeneModel.geneStart and stop > curGeneModel.geneStop) or
				  (start < curGeneModel.geneStart and stop <= curGeneModel.geneStop)):
				return i
#				return curGeneModel, i
			else:
				# interval after the last gene
				if start > curGeneModel.geneStop:
					createNewGenes(geneModels, len(geneModels)-1, ctg, start, stop)
					return i
#					return curGeneModel, -1
				# intergenic region before the last gene
				elif stop <= curGeneModel.geneStart and start >= preGeneModel.geneStop:
					createNewGenes(geneModels, len(geneModels)-2, ctg, start, stop)
					return i-1
#					return curGeneModel, i-1
		else:
			# interval hits the gene
			if start >= curGeneModel.geneStart and stop <= curGeneModel.geneStop:
				return i
#				return curGeneModel, i
			# interval spans the gene
			elif ((start < curGeneModel.geneStart and stop > curGeneModel.geneStart) or
				  (start >= curGeneModel.geneStart and stop > curGeneModel.geneStop) or
				  (start < curGeneModel.geneStart and stop <= curGeneModel.geneStop)):
				return i
#				return curGeneModel, i
			# interval in between two genes
			else:
				if stop <= curGeneModel.geneStart and start >= preGeneModel.geneStop:
					return i
#					return curGeneModel, i
		preGeneModel = curGeneModel
#	all of the above conditions should include every cases,
#	alway return this if there is a bug
	return -2			# -2 encodes cases that failed to consider
#	return None, -2

def cleanContigPairs(dContigPairs, dGFFs):
	distalPairsfile = "trial2.pair"
	dGenePairs = collections.defaultdict()
	daddedModels = collections.defaultdict(list)
	fOUT = open(distalPairsfile, 'w')
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
			geneModelA, geneIndexA = None, -1
			geneModelB, geneIndexB = None, -1
			if ctgA in dGFFs:
				geneIndexA = matchGene(dGFFs[ctgA], startA, stopA)
#				geneModelA, geneIndexA = matchGene(dGFFs[ctgA], startA, stopA)
			if ctgB in dGFFs:
				geneIndexB = matchGene(dGFFs[ctgB], startB, stopB)
#				geneModelB, geneIndexB = matchGene(dGFFs[ctgB], startB, stopB)

			nGeneToLeftA, nGeneToRightA = -1, -1
			nGeneToLeftB, nGeneToRightB = -1, -1
			# use geneIndex as check conditions
			if geneIndexA != -1 and geneIndexB != -1:
				nGeneToLeftA = geneIndexA
				nGeneToRightA = len(dGFFs[ctgA]) - 1 - geneIndexA
				nGeneToLeftB = geneIndexB
				nGeneToRightB = len(dGFFs[ctgB]) - 1 - geneIndexB
				# case where RR with 55
				if (senseA == senseB and senseA == '-' and
					nGeneToLeftA <= nGeneToRightA and
					nGeneToLeftB <= nGeneToRightB and
					nGeneToLeftA + nGeneToLeftB <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, geneIndexA, len(dGFFs[ctgA]), geneIndexB, len(dGFFs[ctgB]), nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == senseB and senseA == '+' and
					nGeneToLeftA >= nGeneToRightA and
					nGeneToLeftB >= nGeneToRightB and
					nGeneToRightA + nGeneToRightB <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, geneIndexA, len(dGFFs[ctgA]), geneIndexB, len(dGFFs[ctgB]), nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == '+' and senseB == '-' and
					nGeneToLeftA >= nGeneToRightA and
					nGeneToLeftB <= nGeneToRightB and
					nGeneToRightA + nGeneToLeftB <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == '-' and senseB == '+' and
					nGeneToLeftA <= nGeneToRightA and
					nGeneToLeftB >= nGeneToRightB and
					nGeneToRightB + nGeneToLeftA <= 1):
					print "\t", "both", pairInfo[i]
					print "\t", "both", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				else:
					print "\t", "both else", pairInfo[i], geneIndexA, geneIndexB, nGeneToLeftA, nGeneToRightA, nGeneToLeftB, nGeneToRightB, senseA, senseB
					print "\t", "both else", nGeneToLeftA + nGeneToRightB, nGeneToLeftA <= nGeneToRightA, nGeneToLeftB >= nGeneToRightB
					pairToRemove.append(i)
			elif geneIndexA == -1 and geneIndexB != -1:
				nGeneToLeftB = geneIndexB
				nGeneToRightB = len(dGFFs[ctgB]) - 1 - geneIndexB
				if (senseA == '+' and senseB == '-' and
					nGeneToLeftB <= nGeneToRightB and
					nGeneToLeftB <= 1):
					print "\t", "EitherA", pairInfo[i]
					print "\t", "EitherA", ctgA, ctgB, "None", geneModelB.geneID, "N/A", "N/A", nGeneToLeftB, nGeneToRightB, senseA, senseB
#					geneModelA, daddedModels = addGeneModel(ctgA, daddedModels, startA, stopA)
#					print geneModelA.lcds, daddedModels[ctgA][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == '-' and senseB == '+' and
					  nGeneToLeftB >= nGeneToRightB and
					  nGeneToRightB <= 1):
					print "\t", "EitherA", pairInfo[i]
					print "\t", "EitherA", ctgA, ctgB, "None", geneModelB.geneID, "N/A", "N/A", nGeneToLeftB, nGeneToRightB, senseA, senseB
#					geneModelA, daddedModels = addGeneModel(ctgA, daddedModels, startA, stopA)
#					print geneModelA.lcds, daddedModels[ctgA][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == senseB and senseA == '+' and
					  nGeneToLeftB > nGeneToRightB and
					  nGeneToRightB <= 1):
					print "\t", "EitherA", pairInfo[i]
					print "\t", "EitherA", ctgA, ctgB, "None", geneModelB.geneID, "N/A", "N/A", nGeneToLeftB, nGeneToRightB, senseA, senseB
#					geneModelA, daddedModels = addGeneModel(ctgA, daddedModels, startA, stopA)
#					print geneModelA.lcds, daddedModels[ctgA][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == senseB and senseA == '-' and
					  nGeneToLeftB < nGeneToRightB and
					  nGeneToLeftB <= 1):
					print "\t", "EitherA", pairInfo[i]
					print "\t", "EitherA", ctgA, ctgB, "None", geneModelB.geneID, "N/A", "N/A", nGeneToLeftB, nGeneToRightB, senseA, senseB
#					geneModelA, daddedModels = addGeneModel(ctgA, daddedModels, startA, stopA)
#					print geneModelA.lcds, daddedModels[ctgA][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				else:
					print "\t", "Either A else", pairInfo[i], geneIndexA, geneIndexB, nGeneToLeftB, nGeneToRightB, senseA, senseB
					print "\t", "Either A else", nGeneToLeftB <= nGeneToRightB, nGeneToRightB
					pairToRemove.append(i)
			elif geneIndexA != -1 and geneIndexB == -1:
				nGeneToLeftA = geneIndexA
				nGeneToRightA = len(dGFFs[ctgA]) - 1 - geneIndexA
				if (senseA == '+' and senseB == '-' and
					nGeneToLeftA >= nGeneToRightA and
					nGeneToRightA <= 1):
					print "\t", "EitherB", pairInfo[i]
					print "\t", "EitherB", ctgA, ctgB, geneModelA.geneID, "None", nGeneToLeftA, nGeneToRightA, "N/A", "N/A", senseA, senseB
#					geneModelB, daddedModels = addGeneModel(ctgB, daddedModels, startB, stopB)
#					print geneModelB.lcds, daddedModels[ctgB][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == '-' and senseB == '+' and
					  nGeneToLeftA <= nGeneToRightA and
					  nGeneToLeftA <= 1):
					print "\t", "EitherB", pairInfo[i]
					print "\t", "EitherB", ctgA, ctgB, geneModelA.geneID, "None", nGeneToLeftA, nGeneToRightA, "N/A", "N/A", senseA, senseB
#					geneModelB, daddedModels = addGeneModel(ctgB, daddedModels, startB, stopB)
#					print geneModelB.lcds, daddedModels[ctgB][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == senseB and senseA == '+' and
					  nGeneToLeftA > nGeneToRightA and
					  nGeneToRightA <= 1):
					print "\t", "EitherB", pairInfo[i]
					print "\t", "EitherB", ctgA, ctgB, geneModelA.geneID, "None", nGeneToLeftA, nGeneToRightA, "N/A", "N/A", senseA, senseB
#					geneModelB, daddedModels = addGeneModel(ctgB, daddedModels, startB, stopB)
#					print geneModelB.lcds, daddedModels[ctgB][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				elif (senseA == senseB and senseA == '-' and
					  nGeneToLeftA < nGeneToRightA and
					  nGeneToLeftA <= 1):
					print "\t", "EitherB", pairInfo[i]
					print "\t", "EitherB", ctgA, ctgB, geneModelA.geneID, "None", nGeneToLeftA, nGeneToRightA, "N/A", "N/A", senseA, senseB
#					geneModelB, daddedModels = addGeneModel(ctgB, daddedModels, startB, stopB)
#					print geneModelB.lcds, daddedModels[ctgB][0].lcds
#					if (geneModelA, geneModelB) not in dGenePairs:
#						dGenePairs[geneModelA, geneModelB] = 1
					fOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(readID, ctgA, startA, senseA, ctgB, startB, senseB))
				else:
					print "\t", "Either B else", pairInfo[i], geneIndexA, geneIndexB, nGeneToLeftA, nGeneToRightA, senseA, senseB
					print "\t", "Either B else", nGeneToLeftA > nGeneToRightA, nGeneToRightA
					pairToRemove.append(i)
			elif geneIndexA == -1 and geneIndexB == -1:
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
			print len(pairToRemove)
			print len(pairInfo)
			# here I need to update the daddedModels
			# to get rid of overlapping exons
			tmp = []
			for i in xrange(len(pairInfo)):
				if i not in pairToRemove:
					tmp.append(pairInfo[i])
				else:
					pairToRemove.remove(i)
			dContigPairs[ctgPair] = tmp
	fOUT.close()
	return distalPairsfile

def rnaPath(icontig, isam, igff, outpathfilename,
			outcontigfilename, pathPrefix="RNAPATH", overlap=30):

	outpathfile = open(outpathfilename, "w")

	outheader = "#settings: %s" % " ".join(sys.argv)
	print outheader
	print >> outpathfile, outheader

	contigNum, nameList, contigDict, origSize = getContigsFromFile(icontig)
	halfSize = calculateN50(origSize)
	print "Getting Gene models ..."
	dGFFs = gff2GeneModels(igff)
	print "Getting discordant contig pairs from SAM ..."
	dContigPairs = asam.sam2CtgPairs(isam, 5)

	distalPairsfile = cleanContigPairs(dContigPairs, dGFFs)
	sys.exit()
	print "building the adjacency graph"

	pathList, edgeSenseDict, visitedDict = getPath(contigNum, distalPairsfile, nameList)

	print pathList
	print "found %d paths" % len(pathList)

	newSizeList = []
	pathID = 0
	outcontigfile = open(outcontigfilename, "w")
	for path in pathList:
		pathID += 1
		# pathPrefix by default "RNAPATH"
		outpathfile.write("%s%d: %s\n" % (pathPrefix, pathID, str(path)))
		vertexNameList = []
		for vertex in path:
			vertexNameList.append(nameList[vertex])
			pathDescription = string.join(vertexNameList, ",")

		print >> outpathfile, pathDescription
		currentVertex = path[0]
		currentSense = "+"
		assemblyList = currentVertex
		sequence = contigDict[currentVertex]
		for nextVertex in path[1:]:
			if (currentVertex, nextVertex) in edgeSenseDict:
				senseList = edgeSenseDict[currentVertex, nextVertex]
				FR = senseList.count(("+", "-"))
				RF = senseList.count(("-", "+"))
			else:
				senseList = edgeSenseDict[nextVertex, currentVertex]
				# flip because the from and end vertices fliped
				FR = senseList.count(("-", "+"))
				RF = senseList.count(("+", "-"))

			FF = senseList.count(("+", "+"))
			RR = senseList.count(("-", "-"))

			if currentSense == "-":
				# we had flipped the upstream piece! Must flip again
				temp1 = FR
				temp2 = FF
				FR = RR
				FF = RF
				RR = temp1
				RF = temp2

			if FR >= FF and FR >= RR and FR >= RF:
				# we have FR - leave alone
				sense1 = "+"
				sense2 = "-"
				assemblyList = ((assemblyList, "+"), (nextVertex, "+"))
				seqleft = sequence[-20:]
				# overlap by default: 30
				seqright = contigDict[nextVertex][:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (currentVertex, nextVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence += contigDict[nextVertex][offset:]
				else:
					sequence += "NN" + contigDict[nextVertex]
				currentSense = "+"
			elif FF >= RR and FF >= RF:
				# we have FF - flip seqright
				sense1 = "+"
				sense2 = "+"
				assemblyList = ((assemblyList, "+"), (nextVertex, "-"))
				seqleft = sequence[-20:]
				seqright = complement(contigDict[nextVertex])[:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (nextVertex, currentVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence += complement(contigDict[nextVertex])[offset:]
				else:
					sequence += "NN" + complement(contigDict[nextVertex])
				currentSense = "-"
			elif RR >= RF:
				# we have RR - flip seqleft
				sense1 = "-"
				sense2 = "-"
				assemblyList = ((assemblyList, "-"), (nextVertex, "+"))
				seqleft = complement(sequence)[:20]
				seqright = contigDict[nextVertex][:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (nextVertex, currentVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence = complement(sequence) + contigDict[nextVertex][offset:]
				else:
					sequence = complement(sequence) + "NN" + contigDict[nextVertex]
				currentSense = "+"
			else:
				# we have RF - flip both
				sense1 = "-"
				sense2 = "+"
				assemblyList = ((assemblyList, "-"), (nextVertex, "-"))
				seqleft = complement(sequence)[-20:]
				seqright = complement(contigDict[nextVertex])[:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (nextVertex, currentVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence = complement(sequence) + complement(contigDict[nextVertex])[offset:]
				else:
					sequence = complement(sequence) + "NN" + complement(contigDict[nextVertex])
				currentSense = "-"

			outstring = "(%d, %d): FF %d RR %d RF %d FR %d : %s %s\t%s" % (currentVertex, nextVertex, FF, RR, RF, FR, sense1, sense2, str(assemblyList))
			print outstring
			print >> outpathfile, outstring
			currentVertex = nextVertex

		outcontigfile.write(">%s%d %dbp %s | %s\n%s\n" % (pathPrefix, pathID, len(sequence), pathDescription, str(assemblyList), sequence))
		newSizeList.append(len(sequence))

	for vertex in contigDict:
		if vertex in visitedDict:
			continue

		newSizeList.append(len(contigDict[vertex]))
		outcontigfile.write(">%s\n%s\n" % (nameList[vertex], contigDict[vertex]))
	calculateN50(newSizeList, referenceMean=halfSize)

def calculateN50(sizeList, referenceMean=None):
	if referenceMean is None:
		totalSize = sum(sizeList)
		referenceMean = totalSize / 2

	sizeList.sort()
	sizeList.reverse()
	currentTotalLength = 0
	for size in sizeList:
		if currentTotalLength + size > referenceMean:
			print "#contigs", len(sizeList)
			print "N50", size
			break

		currentTotalLength += size

#	print sizeList[:50]

	return referenceMean

def getContigsFromFile(contigFileName):
	try:
		incontigfile = open(contigFileName)
	except IOError:
		print "Error opening contig file: %s" % contigFileName

	seq = ""
	contigNum = 0
	nameList = []
	origSize = []
	contigDict = {}
	with open(contigFileName, 'r') as incontigfile:
		for line in incontigfile:
			if line.startswith('>'):
				if seq != "":
					contigDict[contigNum] = seq
					origSize.append(len(seq))
					contigNum += 1
					seq = ""
				chrom = line.strip()[1:]
				nameList.append(chrom)
			else:
				seq += line.strip()

	contigDict[contigNum-1]=seq
	origSize.append(len(seq))
	incontigfile.close()
	return contigNum, nameList, contigDict, origSize

def getPath(contigNum, distalPairsfile, nameList):
	edgeMatrix = EdgeMatrix(contigNum)

	print len(edgeMatrix.edgeArray)
	try:
		print len(edgeMatrix.edgeArray[50])
	except IndexError:
		pass

	print "processing distal pairs"
	verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = processDistalPairsFile(distalPairsfile, edgeMatrix, nameList)

#	print "notSoloDict",notSoloDict #added
#	print "verticesWithEdges",verticesWithEdges #added

	willVisitList = verticesWithEdges.keys()
	willVisitList.sort() #willVisitList before [0, 1, 2, 3, 4, 5, 6]
	print "visiting %d vertices" % len(willVisitList)

	print "cleaning up graph of edges with weight 1"
	verticesToDelete = []
#	print "edgeMatrix before", edgeMatrix.edgeArray #added

#	print "vertexEdges",vertexEdges #added, vertexEdges {0: [3, 4, 5, 6], 1: [2], 2: [1], 3: [0], 4: [0], 5: [0], 6: [0]}

	for rindex in willVisitList: # if a contig is not in notSoloDict, it means that all connections to/from this contig have weight < 1(cutoff), and those weight will be converted to Zero in edgeMatrix.edgeArray
		if rindex not in notSoloDict:
			cindex = vertexEdges[rindex][0]
			edgeMatrix.edgeArray[rindex][cindex] = 0
			edgeMatrix.edgeArray[cindex][rindex] = 0
			verticesToDelete.append(rindex)
#	print "edgeMatrix.edgeArray after 1", edgeMatrix.edgeArray

	for vertex in verticesToDelete:
		willVisitList.remove(vertex)   
	print "%d 1-edges zeroed out" % len(verticesToDelete)

	zeroedEdge = 0
	print "visiting %d vertices" % len(willVisitList)

	leafList = []
#	print "willVisitList after",willVisitList
	print "picking top 2 edges per vertex - zero out others"
	for rindex in willVisitList:
		vertices = vertexEdges[rindex]
		rEdges = []
		for avertex in vertices:
			if avertex in willVisitList:
				rEdges.append((edgeMatrix.edgeArray[rindex][avertex], avertex))

		if len(rEdges) > 2:
			rEdges.sort(reverse=True)
#			rEdges.reverse()
			zeroedEdge += len(rEdges[2:])
			for (weight, cindex) in rEdges[2:]:
				edgeMatrix.edgeArray[rindex][cindex] = 0  #further zero out the non-top2-weight edges
				edgeMatrix.edgeArray[cindex][rindex] = 0
			# this is not necessary
			leafList.append(rindex) #added in RNAPATH*
		elif len(rEdges) == 1:
			if edgeMatrix.edgeArray[rindex][rEdges[0][1]] > 1:
				leafList.append(rindex)
		# I don't think this is right
		else:
			leafList.append(rindex) #added in RNAPATH*

#	print "leafList", leafList #added
#	print "edgeMatrix.edgeArray",edgeMatrix.edgeArray #added

	print "zeroed out %d lower-weight edges at vertices with degree > 2" % zeroedEdge
	pathList, visitedDict = traverseGraph(leafList, edgeMatrix)

	return pathList, edgeSenseDict, visitedDict

def traverseGraph(leafList, edgeMatrix):
	pathList = []
	visitedDict = {}
	leafList.sort()
	print "traveling through the graph"
	for rindex in leafList:
		if visitedDict.has_key(rindex):
			pass
		else:
			path = edgeMatrix.visitLink(rindex) # orig
			# same as original line
#			path = edgeMatrix.visitLink(rindex,visitedDict.keys()) #added
			if len(path) > 1:
				for vertex in path:
					visitedDict[vertex] = ""
				print path
				pathList.append(path)
	return pathList, visitedDict

def processDistalPairsFile(distalPairsfilename, edgeMatrix, nameList):
	contigToRowLookup = {}
	verticesWithEdges = {}
	vertexEdges = {}
	notSoloDict = {}
	edgeSenseDict = {}

	distalPairs = open(distalPairsfilename)
	for line in distalPairs:
		if line[0] == "#":
			continue

		fields = line.strip().split()
		contA = "%s" % fields[1]
		try:
			contig1 = contigToRowLookup[contA]
		except KeyError:
			try:
				contig1 = nameList.index(contA)
				contigToRowLookup[contA] = contig1
			except ValueError:
				print "problem with end1: ", line
				continue

		sense1 = fields[3]

		contB = "%s" % fields[4]
		try:
			contig2 = contigToRowLookup[contB]
		except KeyError:
			try:
				contig2 = nameList.index(contB)
				contigToRowLookup[contB] = contig2
			except ValueError:
				print "problem with end2: ", line
				continue

		sense2 = fields[6]

		edgeMatrix.edgeArray[contig1][contig2] += 1
		edgeMatrix.edgeArray[contig2][contig1] += 1
		verticesWithEdges[contig1] = ""
		verticesWithEdges[contig2] = ""
		if (contig1, contig2) in edgeSenseDict:
			edgeSenseDict[contig1, contig2].append((sense1, sense2))
		elif (contig2, contig1) in edgeSenseDict:
			edgeSenseDict[contig2, contig1].append((sense2, sense1))
		else:
			edgeSenseDict[contig1, contig2] = [(sense1, sense2)]

		if contig1 in vertexEdges:
			if contig2 not in vertexEdges[contig1]:
				vertexEdges[contig1].append(contig2)
		else:
			vertexEdges[contig1] = [contig2]

		if contig2 in vertexEdges:
			if contig1 not in vertexEdges[contig2]:
				vertexEdges[contig2].append(contig1)
		else:
			vertexEdges[contig2] = [contig1]

		if edgeMatrix.edgeArray[contig1][contig2] > 1:
			notSoloDict[contig1] = "" #if a contig has at least one connections with weight >1, it will be added to the notSoloDict
									  # in other words, if all connections to/from a contig have weight <1, it will not be visited in scaffolding path
			notSoloDict[contig2] = ""

	distalPairs.close()

	return verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict

class EdgeMatrix:
	""" Describes a sparse matrix to hold edge data.
	"""

	def __init__(self, dimension):
		self.dimension = dimension
		self.edgeArray = zeros((self.dimension, self.dimension), int16)

	def visitLink(self, fromVertex, ignoreList=[]):
		returnPath = [fromVertex]
		toVertex = []
		for toindex in xrange(self.dimension):
			if self.edgeArray[fromVertex][toindex] > 1 and toindex not in ignoreList:
				toVertex.append(toindex)

		for vertex in toVertex:
			# this is the base case where no further extension can be found
			if sum(self.edgeArray[vertex]) == self.edgeArray[fromVertex][vertex]:
				self.edgeArray[fromVertex][vertex] = 0
				self.edgeArray[vertex][fromVertex] = 0
				return returnPath + [vertex]
			else:
				self.edgeArray[fromVertex][vertex] = 0
				try:
					return returnPath + self.visitLink(vertex, returnPath)
				except IOError:
					return returnPath + [vertex]
		return []

if __name__ == "__main__":
    main(sys.argv)
