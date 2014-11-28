#!/usr/bin/python
import sys

# sample input:scaffold_14442  38932	38983	scaffold_9680	1154	1220	SRR350959.18	12	-	-
# sample ouput:scaffold_14442  38983	38932	scaffold_9680	1220	1154

while True:
	line = sys.stdin.readline()
	if not line:
		break
	lineList = line.strip("\n").split()
	ctg1 = lineList[0]
	ctg2 = lineList[3]
	sense1 = lineList[-2]
	sense2 = lineList[-1]
	if "-" in sense1:
		start1 = lineList[2]
		end1 = lineList[1]
	else:
		start1 = lineList[1]
		end1 = lineList[2]
	if "-" in sense2:
		start2 = lineList[5]
		end2 = lineList[4]
	else:
		start2 = lineList[4]
		end2 = lineList[5]
	outline = "%s\t%s\t%s\t%s\t%s\t%s"%(ctg1,start1,end1,ctg2,start2,end2)
	print outline

