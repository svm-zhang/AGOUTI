#!/usr/bin/python
"""
SSPACE is designed for scaffolding with DNAseq data and GAP size is estimated
based on insert size of PE reads. When using RNAseq data for scaffolding, 
the GAP size would be incorrectly calculated by SSPACE.
This script replaced the incorrect GAP (long stretches of 'N's in output fasta 
with as fixed length of gap containing 100 'N's

last-modified: 10-9-2014
"""

import os
import string
import re
import sys
import optparse


use_message = '''
Replace gaps created by SSPACE into a string of 100 "N's as the putative gap size.

Usage:
    python %prog  <sspaceout.scaffold.fasta>  <outDir> 

'''

def main(argv=None):
    if not argv:
        argv = sys.argv 
    infileDir = argv[1]
    outDir = argv[2]
    path,infileName = os.path.split(infileDir)
    outName = "gapparsed-" + infileName
    outfileName = os.path.join(outDir,outName)

    contigNum, nameList, contigDict_before,nameDict = loadFASTA(infileDir) 
    print "input FASTA file loaded, now converting gaps..."
    contigDict_after = replaceWrongGap(contigDict_before)
    print "gaps converted, now writing output FASTA file..."
    writeOutFile(outfileName,contigDict_after,nameList,nameDict)
    print "output FASTA file created"

def loadFASTA(contigFileName):
    nameList = []
    nameDict = {}
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
                
            #print chrom #added
            nameList.append(chrom)
            nameDict[chrom] = line
            contigNum += 1
            contigDict[chrom] = ""
            if seq :
                prevChrom = nameList[contigNum-2]
                contigDict[prevChrom]=seq
                #origSize.append(len(seq))
                seq=""
        else:
            seq += line.strip()

    contigDict[chrom]=seq
    #origSize.append(len(seq))
    incontigfile.close()
    return contigNum, nameList, contigDict, nameDict

def replaceWrongGap(contigDict):
    """
    Replace gaps created by SSPACE into a string of 'N's with fixed size.
    """
    pattern_old = re.compile("N{1,}")
    pattern_new = "N"*100
    for key in contigDict:
        oldSeq = contigDict[key]
        newSeq = re.sub(pattern_old,pattern_new,oldSeq)
        contigDict[key] = newSeq
    return contigDict

def writeOutFile(outfileName, contigDict_after,nameList,nameDict):
    """
    """
    outfile = open(outfileName,"w")
    for name in nameList:
        outlineName = outfile.write(nameDict[name])
        outlineSeq = outfile.write("%s\n"%(contigDict_after[name]))
    outfile.close()    


if __name__ == "__main__":
    main()
