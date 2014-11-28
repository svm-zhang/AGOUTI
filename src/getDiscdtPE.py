#!/usr/bin/python

import sys
import math
import re
import random
import os

use_message = '''

Usage:
        python %prog  <in.q5.aln> <out.discdt.q5.aln>

'''

if len(sys.argv) != 3:
    print use_message
    sys.exit(1)

inAln = sys.argv[1]
outAln = sys.argv[2]


# since reads with mapQ<5 has been removed, some read pairs 
# has only R1 or R2 left; find out these singletons and remove them:
alnCount = {}
for line in open(inAln,"r"):
    readName = line.split()[0]
    try:
        alnCount[readName] += 1
    except KeyError: 
        alnCount[readName] = 0
        alnCount[readName] += 1

outFile = open(outAln,"w")
for line in open(inAln):
    if alnCount[line.split()[0]] == 2 :
        outFile.write(line)
outFile.close()
