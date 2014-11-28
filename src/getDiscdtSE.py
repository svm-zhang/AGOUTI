#!/usr/bin/python

import sys
import math
import re
import random
import os

use_message = '''

Usage:
        python getDiscdtSE.py  <in.discdt.txt> 

'''
if len(sys.argv) != 2:
    print use_message
    sys.exit(1)
inAln = sys.argv[1]

# readthrough the sam lines with"SA:Z"(split alignment) tag and 
# count the alignment of each reads
alnCount = {}
for line in open(inAln,"r"):
    readName = line.split()[0]
    if readName not in alnCount:
        alnCount[readName] = 0
    alnCount[readName] += 1

discdtReads = {}

for line in open(inAln):
    splitLine = line.split()
    readName = splitLine[0]
    if alnCount[readName] != 2:
        continue
    if readName not in discdtReads:
        discdtReads[readName] = set()
    discdtReads[readName].add(splitLine[2])

keptRead = [key for key in discdtReads if len(discdtReads[key]) == 2]
keptRead = set(keptRead)

for line in open(inAln):
    if (line.split()[0]) in keptRead:
        print line,
