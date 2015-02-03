import os
import sys
import collections
import shlex
import subprocess
#from scipy import stats
#import numpy as np

parpath = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), os.pardir))
if os.path.exists(os.path.join(parpath, "lib", "__init__.py")):
	sys.path.insert(0, parpath)

from lib import agouti_gff as agff

class SZ_ContigPair(object):
	def __init__(self):
		self.ctgA = ""
		self.ctgB = ""
		self.geneA = ""
		self.geneB = ""
		self.exonIndexA = 0
		self.exonIndexB = 0
		self.intervalA = ()
		self.intervalB = ()

	def setContigID(self, ctgA, ctgB):
		self.ctgA = ctgA
		self.ctgB = ctgB

	def setInterval(self, intervalA, intervalB):
		self.intervalA = intervalA
		self.intervalB = intervalB

	def setGene(self, geneA, geneB):
		self.geneA = geneA
		self.geneB = geneB

	def setExonIndex(self, exonIndexA, exonIndexB):
		self.exonIndexA = exonIndexA
		self.exonIndexB = exonIndexB

	def setID(self, geneA, geneB):
		self.geneA = geneA
		self.geneB = geneB

def func_getGFF(iGeneModel):
	dGFFs = collections.defaultdict(list)
	nGene = 0
	with open(iGeneModel, 'r') as fIN:
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
						dGFFs[lobj_GeneModels[geneIndex-1].ctgID].append(lobj_GeneModels[geneIndex-1])
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
	for k, v in dGFFs.iteritems():
		nGeneModels += len(v)
	print "Number of Gene Models parsed: %d" %(nGeneModels)
	return dGFFs

def func_findOverlap(ctgID, start, stop, lobj_GeneModels):
	for i in range(len(lobj_GeneModels)):
		if (start > lobj_GeneModels[i].lcds[-1] or
			lobj_GeneModels[i].lcds[0] > stop):
			continue
		exonIndex = 0
		for j in range(0, len(lobj_GeneModels[i].lcds), 2):
			exonIndex += 1
			if start >= lobj_GeneModels[i].lcds[j] and stop <= lobj_GeneModels[i].lcds[j+1]:
				if not lobj_GeneModels[i].is_fullGene():
					# the value of exonIndex neneeds to be checked
					return lobj_GeneModels[i], exonIndex
	return None, 0

def func_ctgPair2genePair(iContigPair, dGFFs):
	nSuccess = 0
	dGenePairsLinked = collections.defaultdict(list)
	lobj_ctgPairs = []
	with open(iContigPair, 'r') as fIN:
		for line in fIN:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				obj_GeneModelA, exonIndexA = func_findOverlap(tmp_line[0], int(tmp_line[1]), int(tmp_line[2]), dGFFs[tmp_line[0]])
				obj_GeneModelB, exonIndexB = func_findOverlap(tmp_line[3], int(tmp_line[4]), int(tmp_line[5]), dGFFs[tmp_line[3]])
				obj_ContigPair = SZ_ContigPair()
				if (obj_GeneModelA is not None and
				    obj_GeneModelB is not None):
					if exonIndexA and exonIndexB:
						# model A and B must be complementary to each other
						# Both of them cannot have either start or stop codon at the same time
						if ((obj_GeneModelA.startCodon and obj_GeneModelB.stopCodon) or
							(obj_GeneModelA.startCodon and obj_GeneModelB.missStartStop()) or
							(obj_GeneModelA.missStartStop() and obj_GeneModelB.stopCodon) or
							(obj_GeneModelA.stopCodon and obj_GeneModelB.startCodon) or
							(obj_GeneModelA.stopCodon and obj_GeneModelB.missStartStop()) or
							(obj_GeneModelA.missStartStop() and obj_GeneModelB.startCodon) or
							(obj_GeneModelA.missStartStop() and obj_GeneModelB.missStartStop()):
							obj_ContigPair.setContigID(tmp_line[0], tmp_line[3])
							obj_ContigPair.setInterval((tmp_line[1], tmp_line[2]), (tmp_line[4], tmp_line[5]))
							obj_ContigPair.setGene(obj_GeneModelA.geneID, obj_GeneModelB.geneID)
							obj_ContigPair.setExonIndex(exonIndexA, exonIndexB)
							lobj_ctgPairs.append(obj_ContigPair)
							if not obj_GeneModelA.geneID in dGenePairsLinked:
								dGenePairsLinked[obj_GeneModelA.geneID] = [(obj_GeneModelB.geneID, 1, obj_GeneModelA.lcds, obj_GeneModelB.lcds)]
								nSuccess += 1
							else:
								exist = 0
								v = dGenePairsLinked[obj_GeneModelA.geneID]
								for i in range(len(v)):
									if obj_GeneModelB.geneID == v[i][0]:
										v[i] = (obj_GeneModelB.geneID, v[i][1]+1, obj_GeneModelA.lcds, obj_GeneModelB.lcds)
										exist = 1
										break
								if not exist:
									dGenePairsLinked[obj_GeneModelA.geneID] += [(obj_GeneModelB.geneID, 1, obj_GeneModelA.lcds, obj_GeneModelB.lcds)]
#							print line.strip()
	print len(dGenePairsLinked)
	print len(lobj_ctgPairs)
	return dGenePairsLinked, lobj_ctgPairs

def func_filterByLinkNum(dGenePairsLinked, lobj_ctgPairs):
	nLinkFailed = 0
	lobj_filter = []
	print "#####################"
	for i in range(len(lobj_ctgPairs)):
		if lobj_ctgPairs[i].geneA in dGenePairsLinked:
			v = dGenePairsLinked[lobj_ctgPairs[i].geneA]
			for j in range(len(v)):
				if lobj_ctgPairs[i].geneB == v[j][0]:
					print lobj_ctgPairs[i].ctgA, lobj_ctgPairs[i].ctgB, lobj_ctgPairs[i].geneA, lobj_ctgPairs[i].geneB, v[j][1]
					if v[j][1] < 2:
						nLinkFailed += 1
						v.pop(j)
						lobj_filter.append(i)
						break
	for i in range(len(lobj_filter)):
		lobj_ctgPairs.pop(lobj_filter[i])
	print len(dGenePairsLinked)
	print len(lobj_ctgPairs)
	print "Number of exon pairs failed the filter of link number: %d" %(nLinkFailed)
	return dGenePairsLinked, lobj_ctgPairs

def func_prepBED(outp, dGenePairsLinked):
	'''
		Prepare BED file for bedtools coverage
	'''
	lgenes = []
	foGFF = open(outp+"_linkpair.bed", 'w')
	for k, v in dGenePairsLinked.iteritems():
		for i in range(len(v)):
			for j in range(2):
				if not v[i][j].geneID in lgenes:
					exonIndex = 0
					for m in range(0, len(v[i][j].lcds), 2):
						exonIndex += 1
						foGFF.write("%s\t%d\t%d\t%s.Exon%d\n" %(v[i][j].ctgID, v[i][j].lcds[m], v[i][j].lcds[m+1], v[i][j].geneID, exonIndex))
					lgenes.append(v[i][j].geneID)
	foGFF.close()

def func_runCoverageBed(iBAM, outp):
	'''
		Run bedtools coverage
	'''
	try:
		cmd = "bedtools coverage -hist -abam %s -b %s" %(iBAM, outp+"_linkpair.bed")
		p = subprocess.Popen(shlex.split(cmd), stdout = open(outp+"_linkpair.coverage", 'w'))
		p_stderr = p.communicate()[1]
	except (OSError, ValueError) as e:
		sys.stderr.write("%s\n" %(e))
		sys.exit(1)

def func_calculateCovPerExon(iBAM, outp, dGenePairsLinked):
	'''
		Calculate Coverage for Exon Pairs
	'''
	if not os.path.exists(outp+"_linkpair.bed"):
		func_prepBED(outp, dGenePairsLinked)
	if not os.path.exists(outp+"_linkpair.coverage"):
		func_runCoverageBed(iBAM, outp)
	dCovPerExon = collections.defaultdict(list)
	with open(outp+"_linkpair.coverage", 'r') as fCOV:
		for line in fCOV:
			tmp_line = line.strip().split("\t")
			if tmp_line[0] != "all":
				if not tmp_line[3] in dCovPerExon:
#					dCovPerExon[tmp_line[3]] = float(int(tmp_line[4]) * int(tmp_line[5])) / int(tmp_line[6])
					dCovPerExon[tmp_line[3]] = [int(tmp_line[4])] * int(tmp_line[5])
				else:
#					dCovPerExon[tmp_line[3]] += float(int(tmp_line[4]) * int(tmp_line[5])) / int(tmp_line[6])
					dCovPerExon[tmp_line[3]] += [int(tmp_line[4])] * int(tmp_line[5])
	for k, v in sorted(dCovPerExon.iteritems()):
		print k, sum(v)/float(len(v))
	return dCovPerExon

def welch_t_test(mu1, s1, N1, mu2, s2, N2):
	'''
		Perform welch's t test given two samples
	'''
	# Construct arrays to make calculations more succint.
	N_i = np.array([N1, N2])
	dof_i = N_i - 1
	v_i = np.array([s1, s2]) ** 2
	# Calculate t-stat, degrees of freedom, use scipy to find p-value.
	t = (mu1 - mu2) / np.sqrt(np.sum(v_i / N_i))
	dof = (np.sum(v_i / N_i) ** 2) / np.sum((v_i ** 2) / ((N_i ** 2) * dof_i))
	p = stats.distributions.t.sf(np.abs(t), dof) * 2
	return t, p

def func_filterByCovPerExon(iBAM, outp, dGenePairsLinked, lobj_ctgPairs):
	'''
		Filter ctg pairs based on coverage
	'''
	dCovPerExon = func_calculateCovPerExon(iBAM, outp, dGenePairsLinked)
	for i in range(len(lobj_ctgPairs)):
		geneA = lobj_ctgPairs[i].geneA
		geneB = lobj_ctgPairs[i].geneB
		exonIndexA = lobj_ctgPairs[i].exonIndexA
		exonIndexB = lobj_ctgPairs[i].exonIndexB
		keyA = geneA + '.Exon' + str(exonIndexA)
		keyB = geneB + '.Exon' + str(exonIndexB)
		if keyA in dCovPerExon and keyB in dCovPerExon:
			pval = welch_t_test(np.mean(dCovPerExon[keyA]), np.std(dCovPerExon[keyA]), len(dCovPerExon[keyA]),
								np.mean(dCovPerExon[keyB]), np.std(dCovPerExon[keyB]), len(dCovPerExon[keyB]))[1]
#			pval = stats.ttest_rel(dCovPerExon[keyA], dCovPerExon[keyB])[1]
			print geneA, exonIndexA, geneB, exonIndexB, pval
			del dCovPerExon[keyA]
			del dCovPerExon[keyB]

def func_outputctgPairs(outp, lobj_ctgPairs):
	with open(outp+".pair", 'w') as fOPAIR:
		for i in range(len(lobj_ctgPairs)):
			fOPAIR.write("%s\t%s\t%s\t" %(lobj_ctgPairs[i].ctgA, lobj_ctgPairs[i].intervalA[0], lobj_ctgPairs[i].intervalA[1]))
			fOPAIR.write("%s\t%s\t%s\n" %(lobj_ctgPairs[i].ctgB, lobj_ctgPairs[i].intervalB[0], lobj_ctgPairs[i].intervalB[1]))

	with open(outp+".lib", 'w') as fOLIB:
		fOLIB.write("Lib1 TAB %s 99999 0.9 FR\n" %(os.path.realpath(outp+".pair")))

def func_outputGFFs(outp, dGFFs, lobj_ctgPairs):
	fOKEEP = open(outp+".unchanged.gff", 'w')
	fOFIX = open(outp+".tofix.gff", 'w')

	lfixctgs = []
	lfixgenes = []
	for i in range(len(lobj_ctgPairs)):
		if not lobj_ctgPairs[i].ctgA in lfixctgs:
			lfixctgs.append(lobj_ctgPairs[i].ctgA)
		if not lobj_ctgPairs[i].geneA  in lfixgenes:
			lfixgenes.append(lobj_ctgPairs[i].geneA)
		if not lobj_ctgPairs[i].ctgB in lfixctgs:
			lfixctgs.append(lobj_ctgPairs[i].ctgB)
		if not lobj_ctgPairs[i].geneB  in lfixgenes:
			lfixgenes.append(lobj_ctgPairs[i].geneB)

	for k, v in sorted(dGFFs.iteritems()):
		for i in range(len(v)):
			if v[i].ctgID in lfixctgs:
				if v[i].geneID in lfixgenes:
					fOFIX.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;link\n"
								 %(v[i].ctgID, v[i].program, v[i].geneStart,
								   v[i].geneStop, v[i].strand, v[i].geneID))
				else:
					fOFIX.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s\n"
								 %(v[i].ctgID, v[i].program, v[i].geneStart,
								   v[i].geneStop, v[i].strand, v[i].geneID))
				for j in range(0, len(v[i].lcds), 2):
					fOFIX.write("%s\t%s\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s\n"
							 %(v[i].ctgID, v[i].program, v[i].lcds[j],
							   v[i].lcds[j+1], v[i].strand, v[i].geneID, v[i].geneID))
			else:
				fOKEEP.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s\n"
							 %(v[i].ctgID, v[i].program, v[i].geneStart,
							   v[i].geneStop, v[i].strand, v[i].geneID))
				for j in range(0, len(v[i].lcds), 2):
					fOKEEP.write("%s\t%s\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s\n"
							 %(v[i].ctgID, v[i].program, v[i].lcds[j],
							   v[i].lcds[j+1], v[i].strand, v[i].geneID, v[i].geneID))
	fOFIX.close()
	fOKEEP.close()

def main():
	iGeneModel = sys.argv[1]
	iContigPair = sys.argv[2]
	iBAM = sys.argv[3]
	outp = sys.argv[4]

	dGFFs = func_getGFF(iGeneModel)
	dGenePairsLinked, lobj_ctgPairs = func_ctgPair2genePair(iContigPair, dGFFs)
	dGenePairsLinked, lobj_ctgPairs = func_filterByLinkNum(dGenePairsLinked, lobj_ctgPairs)
	func_outputctgPairs(outp, lobj_ctgPairs)
	func_outputGFFs(outp, dGFFs, lobj_ctgPairs)
#	func_filterByCovPerExon(iBAM, outp, dGenePairsLinked, lobj_ctgPairs)

if __name__ == "__main__":
	main()
