import os
import sys
import collections

class AGOUTI_GFF(object):
	def __init__(self):
		self.geneID = ""
		self.geneStart = 0
		self.geneStop = 0
		self.program = ""
		self.ctgID = ""
		self.stopCodon = 0
		self.startCodon = 0
		self.strand = ''
		self.lcds = []
		self.link = 0

	def setGene(self, geneID, start, stop, link):
		self.geneID = geneID
		self.geneStart = start
		self.geneStop = stop
		self.link = link

	def setProgram(self, program):
		self.program = program

	def setContigID(self, ctgID):
		self.ctgID = ctgID

	def setStopCodon(self):
		self.stopCodon = 1

	def setStartCodon(self):
		self.startCodon = 1

	def missStartStop(self):
		if not self.startCodon and not self.stopCodon:
			return True

	def setStrand(self, strand):
		self.strand = strand

	def updateCDS(self, cds_start, cds_stop):
		self.lcds += [cds_start, cds_stop]

	def debug(self):
		print self.geneID, self.geneStart, self.geneStop
		print self.ctgID
		print self.strand

	def getNumExons(self):
		return len(self.lcds)/2

	def is_fullGene(self):
		if self.startCodon and self.stopCodon:
			return True
		else:
			return False

def get_gene_models(gff):
	sys.stderr.write("Getting gene models ... ")
	dGFFs = collections.defaultdict(list)
	nGene = 0
	with open(gff, 'r') as fIN:
		for line in fIN:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				if tmp_line[2] == "gene":
					nGene += 1
		lobj_GeneModels = [AGOUTI_GFF() for i in xrange(nGene)]
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
	sys.stderr.write("%d Gene Models parsed\n" %(nGeneModels))
	return dGFFs
