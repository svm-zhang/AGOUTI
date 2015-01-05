import os
import sys
import collections

from agouti_gff import AGOUTI_GFF

def func_getGeneModelstoFix(iGFFtoFix):
	dGFFs = collections.defaultdict(list)
	nGene = 0
	with open(iGFFtoFix, 'r') as fGFF:
		for line in fGFF:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				if tmp_line[2] == "gene":
					nGene += 1
		lobj_GeneModels = [AGOUTI_GFF() for i in xrange(nGene)]
		geneIndex = -1
		link = 0
		fGFF.seek(0)
		for line in fGFF:
			tmp_line = line.strip().split('\t')
			if not line.startswith('#'):
				if tmp_line[2] == "gene":
					geneIndex += 1
					if len(tmp_line[8].split(';')) == 2:
						link = 1
					else:
						link = 0
					if geneIndex == 0:
						lobj_GeneModels[geneIndex].setGene(tmp_line[8].split('=')[1],
														   int(tmp_line[3]),
														   int(tmp_line[4]), link)
					else:
						dGFFs[lobj_GeneModels[geneIndex-1].ctgID].append(lobj_GeneModels[geneIndex-1])
						lobj_GeneModels[geneIndex].setGene(tmp_line[8].split('=')[1],
														   int(tmp_line[3]),
														   int(tmp_line[4]), link)
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
	return dGFFs

def func_getScaffoldsPath(iEvidence, dGFFs, oGFF):
	nGene = 0
	for k, v in dGFFs.iteritems():
		nGene += len(v)
	print "Number of Gene Models parsed: %d" %(nGene)

	fOFIXGFF = open(oGFF, 'w')
	bool_findPath = 0
	first = 0
	shift = 0
	link = 0
	ngeneAGOUTI = 0
	with open(iEvidence, 'r') as fEVIDENCE:
		for line in fEVIDENCE:
			tmp_line = line.strip().split('|')
			if line.startswith('>'):
				scafID = tmp_line[0][1:9] + '_' +tmp_line[0][9:]
				nCtgs = int(tmp_line[-1][4:])
				if nCtgs > 1:
					bool_findPath = 1
			elif len(tmp_line) == 1:		# skip empty line
				continue
			else:
				if bool_findPath:
					ctgStrand = tmp_line[0].split('_')[0]
					ctgID = "scaffold_" + str(int(tmp_line[0].split('_')[1][3:])-1)
					ctgLen = int(tmp_line[1][4:])
					gapLen = 0
					if len(tmp_line) == 4:
						gapLen = int(tmp_line[-1][4:])
					if nCtgs >= 1:
						print ctgID, ctgLen, ctgStrand, gapLen, shift, nCtgs
					if first == 0:
						shift = 0
						first = 1
					if ctgID in dGFFs:
						if ctgStrand == 'r':
							for i in range(len(dGFFs[ctgID])-1, -1, -1):
								print dGFFs[ctgID][i].geneID
								if dGFFs[ctgID][i].strand == '-':
									update_strand = '+'
								else:
									update_strand = '-'
								ori_geneStart = dGFFs[ctgID][i].geneStop
								ori_geneStop = dGFFs[ctgID][i].geneStart

								# this if/else should be merged as:
								# update_geneStart = shift + (ctgLen-ori_geneStart+1)
								# update_geneStop = shift + (ctgLen-ori_geneStop+1)
								if shift == 0:
									update_geneStart = ctgLen-ori_geneStart+1
									update_geneStop = ctgLen-ori_geneStop+1
								else:
									update_geneStart = shift + (ctgLen-ori_geneStart+1)
									update_geneStop = shift + (ctgLen-ori_geneStop+1)

								# difference between in link or not is just the geneID
								# so I think I should be able to fix this
								if not dGFFs[ctgID][i].link:
									link = 0
									print ("%s\tAGOUTI\tgene\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
										   %(scafID, update_geneStart, update_geneStop, update_strand,
										   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
									for j in range(len(dGFFs[ctgID][i].lcds)-1, -2, -2):
										ori_CDSstart = dGFFs[ctgID][i].lcds[j]
										ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
										if shift == 0:
											update_CDSstart = ctgLen-ori_CDSstart+1
											update_CDSstop = ctgLen-ori_CDSstop+1
										else:
											update_CDSstart = shift + (ctgLen-ori_CDSstart+1)
											update_CDSstop = shift + (ctgLen-ori_CDSstop+1)
										print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
											   %(scafID, update_CDSstart, update_CDSstop, update_strand,
											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
								else:
									# merge gene model
									if link == 0:
										ngeneAGOUTI += 1
										dGFFs[ctgID][i].geneID = "agouti_%d" %(ngeneAGOUTI)
										print ("%s\tAGOUTI\tgene\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
											   %(scafID, update_geneStart, update_geneStop, update_strand,
											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
										link = 1
									else:
										dGFFs[ctgID][i].geneID = "agouti_%d" %(ngeneAGOUTI)
									for j in range(len(dGFFs[ctgID][i].lcds)-1, -2, -2):
										ori_CDSstart = dGFFs[ctgID][i].lcds[j]
										ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
										if shift == 0:
											update_CDSstart = ctgLen-ori_CDSstart+1
											update_CDSstop = ctgLen-ori_CDSstop+1
										else:
											update_CDSstart = shift + (ctgLen-ori_CDSstart+1)
											update_CDSstop = shift + (ctgLen-ori_CDSstop+1)
										print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
											   %(scafID, update_CDSstart, update_CDSstop, update_strand,
											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#									if link == 0:
#										link = 1
#										dGFFs[ctgID][i].geneID = "need to fix"
#										print ("%s\tAGOUTI\tgene\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
#											   %(scafID, update_geneStart, update_geneStop, update_strand,
#											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#										for j in range(len(dGFFs[ctgID][i].lcds)-1, -2, -2):
#											ori_CDSstart = dGFFs[ctgID][i].lcds[j]
#											ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
#											if shift == 0:
#												update_CDSstart = ctgLen-ori_CDSstart+1
#												update_CDSstop = ctgLen-ori_CDSstop+1
#											else:
#												update_CDSstart = shift + (ctgLen-ori_CDSstart+1)
#												update_CDSstop = shift + (ctgLen-ori_CDSstop+1)
#											print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
#												   %(scafID, update_CDSstart, update_CDSstop, update_strand,
#												   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#									else:
#										dGFFs[ctgID][i].geneID = "need to fix"
#										for j in range(len(dGFFs[ctgID][i].lcds)-1, -2, -2):
#											ori_CDSstart = dGFFs[ctgID][i].lcds[j]
#											ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
#											if shift == 0:
#												update_CDSstart = ctgLen-ori_CDSstart+1
#												update_CDSstop = ctgLen-ori_CDSstop+1
#											else:
#												update_CDSstart = shift + (ctgLen-ori_CDSstart+1)
#												update_CDSstop = shift + (ctgLen-ori_CDSstop+1)
#											print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
#												   %(scafID, update_CDSstart, update_CDSstop, update_strand,
#												   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#										link = 0
						elif ctgStrand == 'f':
							for i in range(len(dGFFs[ctgID])):
								if dGFFs[ctgID][i].strand == '-':
									update_strand = '-'
								else:
									update_strand = '+'
								ori_geneStart = dGFFs[ctgID][i].geneStart
								ori_geneStop = dGFFs[ctgID][i].geneStop

								# this if/else should be merged as:
								# update_geneStart = shift+ori_geneStart
								# update_geneStop = shift+ori_geneStop
								if shift == 0:
									update_geneStart = ori_geneStart
									update_geneStop = ori_geneStop
								else:
									update_geneStart = shift + ori_geneStart
									update_geneStop = shift + ori_geneStop

								if not dGFFs[ctgID][i].link:
									link = 0
									print ("%s\tAGOUTI\tgene\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
										   %(scafID, update_geneStart, update_geneStop, update_strand,
										   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
									for j in range(0, len(dGFFs[ctgID][i].lcds), 2):
										ori_CDSstart = dGFFs[ctgID][i].lcds[j]
										ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
										if shift == 0:
											update_CDSstart = ori_CDSstart
											update_CDSstop = ori_CDSstop
										else:
											update_CDSstart = shift + ori_CDSstart
											update_CDSstop = shift + ori_CDSstop
										print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
											   %(scafID, update_CDSstart, update_CDSstop, update_strand,
											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
								else:
									# merge gene models
									if link == 0:
										ngeneAGOUTI += 1
										dGFFs[ctgID][i].geneID = "agouti_%d" %(ngeneAGOUTI)
										print ("%s\tAGOUTI\tgene\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
											   %(scafID, update_geneStart, update_geneStop, update_strand,
											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
										link = 1
									else:
										dGFFs[ctgID][i].geneID = "agouti_%d" %(ngeneAGOUTI)
									for j in range(0, len(dGFFs[ctgID][i].lcds), 2):
										ori_CDSstart = dGFFs[ctgID][i].lcds[j]
										ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
										if shift == 0:
											update_CDSstart = ori_CDSstart
											update_CDSstop = ori_CDSstop
										else:
											update_CDSstart = shift + ori_CDSstart
											update_CDSstop = shift + ori_CDSstop
										print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
											   %(scafID, update_CDSstart, update_CDSstop, update_strand,
											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#									if link == 1:
#										dGFFs[ctgID][i].geneID = "need to fix"
#										for j in range(0, len(dGFFs[ctgID][i].lcds), 2):
#											ori_CDSstart = dGFFs[ctgID][i].lcds[j]
#											ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
#											if shift == 0:
#												update_CDSstart = ori_CDSstart
#												update_CDSstop = ori_CDSstop
#											else:
#												update_CDSstart = shift + ori_CDSstart
#												update_CDSstop = shift + ori_CDSstop
#											print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents%s"
#												   %(scafID, update_CDSstart, update_CDSstop, update_strand,
#												   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#										link = 0
#									else:
#										link = 1
#										dGFFs[ctgID][i].geneID = "need to fix"
#										print ("%s\tAGOUTI\tgene\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents=%s"
#											   %(scafID, update_geneStart, update_geneStop, update_strand,
#											   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
#										for j in range(0, len(dGFFs[ctgID][i].lcds), 2):
#											ori_CDSstart = dGFFs[ctgID][i].lcds[j]
#											ori_CDSstop = dGFFs[ctgID][i].lcds[j-1]
#											if shift == 0:
#												update_CDSstart = ori_CDSstart
#												update_CDSstop = ori_CDSstop
#											else:
#												update_CDSstart = shift + ori_CDSstart
#												update_CDSstop = shift + ori_CDSstop
#											print ("%s\tAGOUTI\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds;Parents%s"
#												   %(scafID, update_CDSstart, update_CDSstop, update_strand,
#												   dGFFs[ctgID][i].geneID, dGFFs[ctgID][i].geneID))
					nCtgs -= 1
					shift += ctgLen + gapLen
					if nCtgs == 0:
						bool_findPath = 0
						first = 1
						shift = 0
	fOFIXGFF.close()

def main():
	iGFFtoFix = sys.argv[1]
	iEvidence = sys.argv[2]		# *.final.evidence file from SSPACE
	oGFF = sys.argv[3]

	dGFFs = func_getGeneModelstoFix(iGFFtoFix)
	func_getScaffoldsPath(iEvidence, dGFFs, oGFF)

if __name__ == "__main__":
	main()
