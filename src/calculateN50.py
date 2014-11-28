#!/usr/bin/python

import os
import string
import re
import sys

use_message = '''
Compute the N50 value of a multi-fasta supercontig/scafffold file.

Usage:
	python calculateN50.py  <in.fastafile> 

'''
def main():
	if len(sys.argv) != 2:
		print use_message
		sys.exit()
	infileName = sys.argv[1]

	print "%s loaded..."%(infileName)

	contigNum, nameList, contigDict, origSize = getContigsFromFile(infileName)	
	N50 = calculateN50(origSize,referenceMean=None)
	print "contigNum: ",contigNum
	print "N50: ",N50

def getContigsFromFile(contigFileName):
	nameList = []
	origSize = []
	contigNum = 0 
	contigDict = {}
	seq = ""

	try:
		incontigfile = open(contigFileName)
	except IOError:
		print "Error opening contig file: %s" % contigFileName
		#return contigNum, nameList, contigDict, origSize

	for line in incontigfile:
		if ">" in line:
			if len(line.split()) >1:
				chrom = line.split()[0][1:]
			else :
				chrom = line[1:-1]
			nameList.append(chrom)
			contigNum += 1
			contigDict[chrom] = ""
			if seq :
				prevChrom = nameList[contigNum-2]
				contigDict[prevChrom]=seq
				origSize.append(len(seq))
				seq=""
		else:
			seq += line.strip()

	contigDict[chrom]=seq
	origSize.append(len(seq))
	incontigfile.close()
	return contigNum, nameList, contigDict, origSize
	
def calculateN50(sizeList,referenceMean=None):
	"""calculate the N50 value of the give size list of a contig/scaffold set.
	"""
	N50 = 0
	if referenceMean is None:
		totalSize = sum(sizeList)
		referenceMean = totalSize / 2
	sizeList.sort()
	sizeList.reverse()
	currentTotalLength = 0
	for size in sizeList:
		if currentTotalLength + size > referenceMean:
			break

		currentTotalLength += size
	return size

if __name__ == "__main__":
	main()
