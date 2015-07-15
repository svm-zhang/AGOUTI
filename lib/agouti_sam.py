import os
import sys
import argparse
import re
import collections
import pysam

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

def get_joining_pairs(bamFile, min_nLinks):
	minFracOvl = 0.0
	maxFracMismatch = 1
	minMapQ = 5
#	fSAM = open(isam, 'r')
	sys.stderr.write("Getting joining pairs ... ")
	dContigPairs = collections.defaultdict(list)
	while True:
		pairA = bamFile.readline().strip().split("\t")
		pairB = bamFile.readline().strip().split("\t")
		if len(pairA) == 1 or len(pairB) == 1:
			break
		contigA = pairA[2]
		contigB = pairB[2]
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
#			print pairA
#			print pairB
#			print flagsA
#			print flagsB
#			print nMismatchesA, nMismatchesB
#			print alnLenA, alnLenB
#			print alnLenA/readLenA, alnLenB/readLenB
#			print nMismatchesA/alnLenA, nMismatchesB/alnLenB
#			sys.exit()

			if (min(alnLenA/readLenA, alnLenB/readLenB) >= minFracOvl and				# minimum fraction of overlaps
				max(nMismatchesA/alnLenA, nMismatchesB/alnLenB) <= maxFracMismatch and	# maximum fraction of mismatches
				min(mapQA, mapQB) >= minMapQ):				# minimum mapping quality
				startA = leftMostPosA + 1
				stopA = startA + 1 + int(alnLenA)
				startB = leftMostPosB + 1
				stopB = startB + 1 + int(alnLenB)
#				startA, stopA = getMappedRegionOnContigs(int(pairA[3]), int(alnLenA), flagsA)
#				startB, stopB = getMappedRegionOnContigs(int(pairB[3]), int(alnLenB), flagsB)
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
#				sys.stdout.write("\t".join(pairA)+"\n")
#				sys.stdout.write("\t".join(pairB)+"\n")
#				if contigA <= contigB:
#					sys.stdout.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\n" %(pairA[0], contigA, startA, senseA,
#																  contigB, startB, senseB))
#				else:
#					sys.stdout.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\n" %(pairB[0], contigB, startB, senseB,
#																  contigA, startA, senseA))
#	sys.exit()

	# filter some of the contig pairs who do not
	# have a minimum number of read support
	nCtgPairs = 0
	for k, v in dContigPairs.items():
		if len(v) < min_nLinks:
			del dContigPairs[k]
		else:
			nCtgPairs += len(v)
	sys.stderr.write("%d joining pairs parsed\n" %(nCtgPairs))
	return dContigPairs

#def get_joining_pairs(bamFile, min_nLinks):
#	print bamFile
#	print min_nLinks
#	return get_joining_pairs(bamFile, args.min_nLinks)

#def main():
#	parser = argparse.ArgumentParser()
#	parser.add_argument("-gff", metavar="FILE", dest="igff", required=True, help="gene models in GFF format")
#	parser.add_argument("-mnl", metavar="INT", dest="min_nLinks", default=5, help="minimum number of reads supporting a link between a contig pair")
#	parser.add_argument("sam", nargs='?', help="reads mapping in SAM format")
#	args = parser.parse_args()

#	cleanContigPairs(dContigPairs, dGFFs)

#if __name__ == "__main__":
#	main()
