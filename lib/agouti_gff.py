import os
import collections
import re
import sys

from lib import agouti_log as agLOG

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
		self.fake = 0
		self.merge = 0

	def setGene(self, geneID, start=0, stop=0, fake=0):
		self.geneID = geneID
		self.geneStart = start
		self.geneStop = stop
		self.fake = fake

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
		print("geneID=%s - start=%d - stop=%d - strand=%s - contig=%s"
			  %(self.geneID, self.geneStart, self.geneStop, self.strand, self.ctgID))

	def getNumExons(self):
		return len(self.lcds)/2

	def is_fullGene(self):
		if self.startCodon and self.stopCodon:
			return True
		else:
			return False

def get_gene_models(gff, outDir, prefix, debug=0):
	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_GFFs")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_gff.progressMeter" %(prefix))
	agGFFProgress = agLOG.PROGRESS_METER(moduleName)
	agGFFProgress.add_file_handler(progressLogFile)
	agGFFProgress.logger.info("[BEGIN] Getting gene models")
	dGFFs = collections.defaultdict(list)
	nGene = 0
	with open(gff, 'r') as fIN:
		for line in fIN:
			if line.startswith("##FASTA") or line.startswith("##Fasta"):
					break
			# skip empty lines and lines starting with '#'
			if not line.startswith('#') and len(line.strip()) > 0:
				tmp_line = line.strip().split("\t")
				if tmp_line[2] == "gene":
					nGene += 1
		if nGene == 0:
			agGFFProgress.logger.error("Found zero genes")
			agGFFProgress.logger.error("Please check your GFF file")
			sys.exit(1)
		lobj_GeneModels = [AGOUTI_GFF() for i in range(nGene)]
		geneIndex = -1
		stop = 0
		fIN.seek(0)
		for line in fIN:
			# Stop before getting into FASTA zone
			if line.startswith("##FASTA") or line.startswith("##Fasta"):
					stop = 1
					break
			# skip empty lines and lines starting with '#'
			if not line.startswith('#') and line.strip():
				tmp_line = line.strip().split("\t")
				if tmp_line[2] == "gene":
					geneIndex += 1
					attrs = tmp_line[8].split(';')
					for attr in attrs:
						attrID, attrVal = attr.split('=')
						if attrID == "ID":
							geneID = attrVal
							break
					#m = re.search("(;ID=.+;|ID=.+;|ID=.+|;ID=.+)", tmp_line[8])
					#print m.group()
					#geneID = m.group().strip(';').split('=')[1]
					if geneIndex == 0:
						#lobj_GeneModels[geneIndex].setGene(tmp_line[8].split('=')[1],
						#								   int(tmp_line[3]),
						#								   int(tmp_line[4]))
						lobj_GeneModels[geneIndex].setGene(geneID,
														   int(tmp_line[3]),
														   int(tmp_line[4]))
					else:
						preCtgID = lobj_GeneModels[geneIndex-1].ctgID
						preGeneID = lobj_GeneModels[geneIndex-1].geneID
						dGFFs[preCtgID].append(lobj_GeneModels[geneIndex-1])
						#lobj_GeneModels[geneIndex].setGene(tmp_line[8].split('=')[1],
						#								   int(tmp_line[3]),
						#								   int(tmp_line[4]))
						lobj_GeneModels[geneIndex].setGene(geneID,
														   int(tmp_line[3]),
														   int(tmp_line[4]))
					lobj_GeneModels[geneIndex].setProgram(tmp_line[1])
					lobj_GeneModels[geneIndex].setContigID(tmp_line[0])
					lobj_GeneModels[geneIndex].setStrand(tmp_line[6])
				elif tmp_line[2] == "stop_codon":
					lobj_GeneModels[geneIndex].setStopCodon()
				elif tmp_line[2] == "start_codon":
					lobj_GeneModels[geneIndex].setStartCodon()
				elif tmp_line[2] == "CDS":
					lobj_GeneModels[geneIndex].updateCDS(int(tmp_line[3]), int(tmp_line[4]))
		if not stop and geneIndex >= 0:
			dGFFs[lobj_GeneModels[geneIndex].ctgID].append(lobj_GeneModels[geneIndex])

	if debug:
		debugLogFile = os.path.join(moduleOutDir, "%s.agouti_gff.debug" %(prefix))
		agGFFDebug = agLOG.DEBUG(moduleName, debugLogFile)
		agGFFDebug.debugger.debug("Sequence\tNum_Gene_Models")

	nGeneModels = 0
	for k, v in sorted(dGFFs.items()):
		genes = [(gene.geneStart, gene.geneStop) for gene in v]
		# make sure gene model are in ascending order
		soGenes = sorted(range(len(genes)), key=lambda k:genes[k])
		tmpV = []
		for i in range(len(soGenes)):
			index = soGenes[i]
			tmpV.append(v[index])
		dGFFs[k] = tmpV
		nGeneModels += len(tmpV)
		if debug:
			agGFFDebug.debugger.debug("%s\t%d" %(k, len(tmpV)))

	agGFFProgress.logger.info("%d Gene Models parsed" %(nGeneModels))
	agGFFProgress.logger.info("[DONE]")
	return dGFFs
