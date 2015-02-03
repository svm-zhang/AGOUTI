import os
import sys
import argparse
import re
import collections

parpath = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), os.pardir))
if os.path.exists(os.path.join(parpath, "lib", "__init__.py")):
	sys.path.insert(0, parpath)

from lib import agouti_gff as agff

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

#	if selfStrand != mateStrand:
#		return (paired, proper, selfUnMapped, mateUnMapped,
#				"FR", secondAln, duplicates)
#	else:
#		if selfStrand == 1:
#			return (paired, proper, selfUnMapped, mateUnMapped,
#					"FF", secondAln, duplicates)
#		else:
#			return (paired, proper, selfUnMapped, mateUnMapped,
#					"RR", secondAln, duplicates)

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

def getContigPairs(isam, min_nLinks):
	dContigPairs = collections.defaultdict(list)
	minFracOvl = 0.8
	maxFracMismatch = 0.1
	minMapQ = 5
	fSAM = open(isam, 'r')
	while True:
		pairA = fSAM.readline().strip().split("\t")
		pairB = fSAM.readline().strip().split("\t")
		if len(pairA) == 1 or len(pairB) == 1:
			break
		contigA = pairA[2]
		contigB = pairB[2]
		if pairA[0] == pairB[0] and contigA != contigB:
			alnLenA = getCIGAR(pairA[5])
			alnLenB = getCIGAR(pairB[5])
			readLenA = len(pairA[9])
			readLenB = len(pairB[9])
			nMismatchesA = getMismatches(pairA[11:])
			nMismatchesB = getMismatches(pairB[11:])
			mapQA = int(pairA[4])
			mapQB = int(pairB[4])
			flagsA = explainSAMFlag(int(pairA[1]))
			flagsB = explainSAMFlag(int(pairB[1]))
#			print flagsA
#			print flagsB
#			print nMismatchesA, nMismatchesB
#			print alnLenA, alnLenB
#			print alnLenA/readLenA, alnLenB/readLenB
#			print nMismatchesA/alnLenA, nMismatchesB/alnLenB
			if (min(alnLenA/readLenA, alnLenB/readLenB) >= minFracOvl and				# minimum fraction of overlaps
				max(nMismatchesA/alnLenA, nMismatchesB/alnLenB) <= maxFracMismatch and	# maximum fraction of mismatches
				min(mapQA, mapQB) >= minMapQ):				# minimum mapping quality
				startA = int(pairA[3])
				stopA = startA + int(alnLenA)
				startB = int(pairB[3])
				stopB = startB + int(alnLenB)
#				startA, stopA = getMappedRegionOnContigs(int(pairA[3]), int(alnLenA), flagsA)
#				startB, stopB = getMappedRegionOnContigs(int(pairB[3]), int(alnLenB), flagsB)
				senseA = flagsA[4]
				senseB = flagsB[4]
				if contigA <= contigB:
					if (contigA, contigB) not in dContigPairs:
						dContigPairs[contigA, contigB] = [(startA, startB, stopA, stopB, senseA, senseB, pairA[0])]
					else:
						dContigPairs[contigA, contigB] += [(startA, startB, stopA, stopB, senseA, senseB, pairA[0])]
				else:
					if (contigB, contigA) not in dContigPairs:
						dContigPairs[contigB, contigA] = [(startB, startA, stopB, stopA, senseB, senseA, pairB[0])]
					else:
						dContigPairs[contigB, contigA] += [(startB, startA, stopB, stopA, senseB, senseA, pairB[0])]
#				if contigA <= contigB:
#					sys.stdout.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\n" %(pairA[0], contigA, startA, senseA,
#																  contigB, startB, senseB))
#				else:
#					sys.stdout.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\n" %(pairB[0], contigB, startB, senseB,
#																  contigA, startA, senseA))

	# filter some of the contig pairs who do not
	# have a minimum number of read support
	nCtgPairs = 0
	for k, v in dContigPairs.items():
		if len(v) < min_nLinks:
			del dContigPairs[k]
		else:
			nCtgPairs += len(v)
	sys.stderr.write("Number of Pairs of Contigs: %d\n" %(nCtgPairs))
	return dContigPairs

def getGeneModels(igff):
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

def matchGene(geneModels, start, stop):
	exonIndex = 0
	for i in xrange(len(geneModels)):
		geneModel = geneModels[i]
		# interval hits the gene
		if start >= geneModel.geneStart and stop <= geneModel.geneStop:
			return geneModel
		# interval spans the gene
		elif ((start >= geneModel.geneStart and start < geneModel.geneStop and
			  stop > geneModel.geneStop) or (start < geneModel.geneStart and
			  stop <= geneModel.geneStop and stop > geneModel.geneStart)):
			return geneModel
#		for j in range(0, len(geneModel.lcds), 2):
#			if start >= geneModel.lcds[j] and stop <= geneModel.lcds[j+1]:
#				exonIndex = (j + 2) / 2
#				return geneModel, exonIndex
#			elif (start >= geneModel.lcds[j] and		# cases where the interval
#				  start <= geneModel.lcds[j+1] and		# hits entirely in introns
#				  stop > geneModel.lcds[j+1]):
#				return geneModel, exonIndex
#			elif (start > geneModel.lcds[j+1] and		# cases where the interval
#				  stop < geneModel.lcds[j+2]):			# spans exon/intron junctions
#				return geneModel, exonIndex
	return None

def contigPairToGenePair():
	pass

def cleanContigPairs(dContigPairs, dGFFs):
	for ctgPair, pairInfo in dContigPairs.items():
		ctgA = ctgPair[0]
		ctgB = ctgPair[1]
		keep = []
		print ">", ctgA, ctgB
		for i in xrange(len(pairInfo)):
			startA = pairInfo[i][0]
			stopA = pairInfo[i][2]
			startB = pairInfo[i][1]
			stopB = pairInfo[i][3]
			senseA = pairInfo[i][4]
			senseB = pairInfo[i][5]
			geneModelA = None
			geneModelB = None
			print pairInfo[i]
			if ctgA in dGFFs:
				geneModelA = matchGene(dGFFs[ctgA], startA, stopA)
			if ctgB in dGFFs:
				geneModelB = matchGene(dGFFs[ctgB], startB, stopB)
			if geneModelA is not None and geneModelB is not None:
				nGeneToRightA = len(dGFFs[ctgA]) - dGFFs[ctgA].index(geneModelA)
				nGeneToLeftA = dGFFs[ctgA].index(geneModelA)
				nGeneToLeftB = dGFFs[ctgB].index(geneModelB) - 0
				nGeneToRightB = len(dGFFs[ctgB]) - dGFFs[ctgB].index(geneModelB)
				# the following could be used as how to connect contigs
				if nGeneToLeftA <= nGeneToRightA:
					geneARelativePos = 5
					nGeneSpannedA = nGeneToLeftA
				else:
					geneARelativePos = 3
					nGeneSpannedA = nGeneToRightA
				if nGeneToLeftB <= nGeneToRightB:
					geneBRelativePos = 5
					nGeneSpannedB = nGeneToLeftB
				else:
					geneBRelativePos = 3
					nGeneSpannedB = nGeneToRightB
				print "\t", ctgA, ctgB, geneModelA.geneID, geneModelB.geneID, nGeneSpannedA, geneARelativePos, nGeneSpannedB, geneBRelativePos
		print

	# I need to further filter contig pairs by the number of gene pairs supported by a minimum number of reads
	# I need to figure out the orientation stuff
	# for contig pairs do not hit gene models, I need to create them and find their position






#		sys.exit()
#				if geneModelA is not None:
#					if exonIndexA != 0:
#						print ctgA, pairInfo[i], geneModelA.geneID, "hit exon %d" %(exonIndexA)
#					else:
#						print ctgA, pairInfo[i], geneModelA.geneID, "alternative splicing"
#				else:
#					print "no gene model match found"
#					print ctgA, pairInfo[i]
#			else:
#				# perhaps I need to create new gene models here
#				print "no gene models found for this contig", ctgA
#				sys.exit()
#				if geneModelB is not None:
#					if exonIndexB != 0:
#						print ctgB, pairInfo[i], geneModelB.geneID, "hit exon %d" %(exonIndexB)
#					else:
#						print ctgB, pairInfo[i], geneModelB.geneID, "alternative splicing"
#				else:
#					print "no gene model match found"
#					print ctgB, pairInfo[i]
#			else:
#				# perhaps I need to create new gene models here
#				print "no gene models found for this contig", ctgB
#				sys.exit()

#			if geneModelA is None and geneModelB is None:
#				# cases where I need to create gene models
#				pass
#			else:
#				if (senseA == '+' and senseB == '+' or
#					senseA == '-' and senseB == '-'):
#				if geneModelA is not None and geneModelB is not None:
#					nGeneSpannedA = len(dGFFs[ctgA]) - dGFFs[ctgA].index(geneModelA)
#					nGeneSpannedB = dGFFs[ctgB].index(geneModelB) - 0
#				print nGeneSpannedA, nGeneSpannedB
#					if nGeneSpannedA + nGeneSpannedB <= 1:
#						keep.append(pairInfo[i])
#		if len(keep) > 0:
#			print keep
	pass

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gff", metavar="FILE", dest="igff", required=True, help="gene models in GFF format")
	parser.add_argument("-mnl", metavar="INT", dest="min_nLinks", default=5, help="minimum number of reads supporting a link between a contig pair")
	parser.add_argument("sam", nargs='?', help="reads mapping in SAM format")
	args = parser.parse_args()

	dGFFs = getGeneModels(args.igff)
	dContigPairs = getContigPairs(args.sam, args.min_nLinks)
	cleanContigPairs(dContigPairs, dGFFs)

if __name__ == "__main__":
	main()

