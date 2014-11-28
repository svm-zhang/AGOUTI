#!/usr/bin/python
import pysam
import sys

use_message = '''
Convert the split-mapped single-end .sam file into .bedpe format.

Usage:
	python makeBedpeFromSE.py  <in.sam>  <out.bedpe>

'''

def main():
	if len(sys.argv) != 3:
		print use_message
		sys.exit()
	insamfileName = sys.argv[1]
	outbedpefileName = sys.argv[2]
	
	samfile = pysam.Samfile(insamfileName,"r")
	outbedpe = open(outbedpefileName,"w")

	linkpair_dict = {}
	for aln in samfile:
		readname = aln.qname
		ctg = samfile.getrname(aln.tid)
		start = aln.pos
		end = aln.aend
		sense = getSense(aln.is_reverse)
		pairInfo = [ctg,start,end,sense]
		try:
			linkpair_dict[readname].append(pairInfo)
		except KeyError:
			linkpair_dict[readname] = []
			linkpair_dict[readname].append(pairInfo)
	for pair in linkpair_dict:
		ctg1,start1,end1,sense1 = linkpair_dict[pair][0]
		ctg2,start2,end2,sense2 = linkpair_dict[pair][1]
		
		if sense1 != sense2:
			#print pair,linkpair_dict[pair][0]
			#print pair,linkpair_dict[pair][1]
			continue
		
		#sample line: scaffold_11050  6579    6680    scaffold_11051  51      152     SRR919327.10    60      +       -
		bedpe = "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n"%\
				(ctg1,start1,end1,ctg2,start2,end2,pair,60,"+","-")
		outbedpe.write(bedpe)

	outbedpe.close()

def getSense(status):
	if status:
		return "-"
	if not status:
		return "+"

if __name__== "__main__":
		main()
