import os
import sys

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

