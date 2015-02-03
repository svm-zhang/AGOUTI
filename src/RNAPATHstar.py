#!/usr/bin/python

import sys
import optparse
import string
from numpy import zeros, int16

versionString = "RNAPATH*" # modified from RNAPATH by luting zhuo

def compNT(nt):
	""" returns the complementary basepair to base nt
	"""
	compDict = { "A": "T",
				 "T": "A",
				 "G": "C",
				 "C": "G",
				 "S": "S",
				 "W": "W",
				 "R": "Y",
				 "Y": "R",
				 "M": "K",
				 "K": "M",
				 "H": "D",
				 "D": "H",
				 "B": "V",
				 "V": "B",
				 "N": "N",
				 "a": "t",
				 "t": "a",
				 "g": "c",
				 "c": "g",
				 "n": "n",
				 "z": "z"
	}

	return compDict.get(nt, "N")


def complement(sequence, length=-1):
	""" returns the complement of the sequence.
	"""
	newSeq = ""

	seqLength = len(sequence)

	if length == seqLength or length < 0:
		seqList = list(sequence)
		seqList.reverse()
		return "".join(map(compNT, seqList))

#TODO: this seems to want to deal with case where length is more than
# sequence length except that a negative index on a sequence is fine
# index will only be overrun if length is negative but that case is
# handled above
	for index in range(seqLength - 1,seqLength - length - 1, -1):
		try:
			newSeq += compNT(sequence[index])
		except:
			newSeq += "N"

	return newSeq


def main(argv=None):
	if not argv:
		argv = sys.argv

	usage = "python %prog incontigfile distalPairs outpathfile outcontigfile [--prefix string] [--overlap bp]"

	# parser = getParser(usage)
	# (options, args) = parser.parse_args(argv[1:])

	if len(argv) < 4:
		print usage
		sys.exit(0)

	incontigfilename = argv[1]
	distalPairsfile = argv[2]
	outpathfilename = argv[3]
	outcontigfilename = argv[4]

	rnaPath(incontigfilename, distalPairsfile, outpathfilename,
			outcontigfilename)


# def getParser(usage):
    # parser = optparse.OptionParser(usage=usage)
    # parser.add_option("--prefix", dest="pathPrefix")
    # parser.add_option("--overlap", type="int", dest="overlap")

    # configParser = getConfigParser()
    # section = "RNAPATH"
    # pathPrefix = getConfigOption(configParser, section, "pathPrefix", "RNAPATH")
    # overlap = getConfigIntOption(configParser, section, "overlap", 30)

    # parser.set_defaults(pathPrefix=pathPrefix, overlap=overlap)

    # return parser


def rnaPath(incontigfilename, distalPairsfile, outpathfilename,
            outcontigfilename, pathPrefix="RNAPATH", overlap=30):

	outpathfile = open(outpathfilename, "w")

	outheader = "#settings: %s" % " ".join(sys.argv)
	print outheader
	print >> outpathfile, outheader

	contigNum, nameList, contigDict, origSize = getContigsFromFile(incontigfilename)
	halfSize = calculateN50(origSize)
	print "building the adjacency graph"
	pathList, edgeSenseDict, visitedDict = getPath(contigNum, distalPairsfile, nameList)

	print "found %d paths" % len(pathList)            

	newSizeList = []
	pathID = 0
	outcontigfile = open(outcontigfilename, "w")
	for path in pathList:
		pathID += 1
		outpathfile.write("%s%d: %s\n" % (pathPrefix, pathID, str(path))) 
		vertexNameList = []
		for vertex in path:
			vertexNameList.append(nameList[vertex])
			pathDescription = string.join(vertexNameList, ",")

		print >> outpathfile, pathDescription
		currentVertex = path[0]
		currentSense = "+"
		assemblyList = currentVertex
		sequence = contigDict[currentVertex]
		for nextVertex in path[1:]:
			if (currentVertex, nextVertex) in edgeSenseDict:
				senseList = edgeSenseDict[currentVertex, nextVertex]
				FR = senseList.count(("+", "-"))
				RF = senseList.count(("-", "+"))
			else:
				senseList = edgeSenseDict[nextVertex, currentVertex]
				# flip
				FR = senseList.count(("-", "+"))
				RF = senseList.count(("+", "-"))

			FF = senseList.count(("+", "+"))
			RR = senseList.count(("-", "-"))
			if currentSense == "-":
				# we had flipped the upstream piece! Must flip again
				temp1 = FR
				temp2 = FF
				FR = RR
				FF = RF
				RR = temp1
				RF = temp2

			if FR >= FF and FR >= RR and FR >= RF:
				# we have FR - leave alone
				sense1 = "+"
				sense2 = "-"
				assemblyList = ((assemblyList, "+"), (nextVertex, "+"))
				seqleft = sequence[-20:]
				seqright = contigDict[nextVertex][:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (currentVertex, nextVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence += contigDict[nextVertex][offset:]
				else:
					sequence += "NN" + contigDict[nextVertex]

				currentSense = "+"
			elif FF >= RR and FF >= RF:
				# we have FF - flip seqright
				sense1 = "+"
				sense2 = "+"
				assemblyList = ((assemblyList, "+"), (nextVertex, "-"))
				seqleft = sequence[-20:]
				seqright = complement(contigDict[nextVertex])[:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (nextVertex, currentVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence += complement(contigDict[nextVertex])[offset:]
				else:
					sequence += "NN" + complement(contigDict[nextVertex])

				currentSense = "-"
			elif RR >= RF:
				# we have RR - flip seqleft
				sense1 = "-"
				sense2 = "-"
				assemblyList = ((assemblyList, "-"), (nextVertex, "+"))
				seqleft = complement(sequence)[:20]
				seqright = contigDict[nextVertex][:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (nextVertex, currentVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence = complement(sequence) + contigDict[nextVertex][offset:]
				else:
					sequence = complement(sequence) + "NN" + contigDict[nextVertex]

				currentSense = "+"
			else:
				# we have RF - flip both
				sense1 = "-"
				sense2 = "+"
				assemblyList = ((assemblyList, "-"), (nextVertex, "-"))
				seqleft = complement(sequence)[-20:]
				seqright = complement(contigDict[nextVertex])[:overlap]
				if seqleft in seqright:
					pos = seqright.index(seqleft)
					offset = pos + 20
					outstring = "stitching %d and %d using %d overlap" % (nextVertex, currentVertex, offset)
					print outstring
					print >> outpathfile, outstring
					sequence = complement(sequence) + complement(contigDict[nextVertex])[offset:]
				else:
					sequence = complement(sequence) + "NN" + complement(contigDict[nextVertex])

				currentSense = "-"

			outstring = "(%d, %d): FF %d RR %d RF %d FR %d : %s %s\t%s" % (currentVertex, nextVertex, FF, RR, RF, FR, sense1, sense2, str(assemblyList))
			print outstring
			print >> outpathfile, outstring
			currentVertex = nextVertex

		outcontigfile.write(">%s%d %dbp %s | %s\n%s\n" % (pathPrefix, pathID, len(sequence), pathDescription, str(assemblyList), sequence))
		newSizeList.append(len(sequence))

	for vertex in contigDict:
		if vertex in visitedDict:
			continue

		newSizeList.append(len(contigDict[vertex]))
		outcontigfile.write(">%s\n%s\n" % (nameList[vertex], contigDict[vertex]))

	calculateN50(newSizeList, referenceMean=halfSize)

def calculateN50(sizeList, referenceMean=None):
	if referenceMean is None:
		totalSize = sum(sizeList)
		referenceMean = totalSize / 2

	sizeList.sort()
	sizeList.reverse()
	currentTotalLength = 0
	for size in sizeList:
		if currentTotalLength + size > referenceMean:
			print "#contigs", len(sizeList)
			print "N50", size
			break

		currentTotalLength += size

	print sizeList[:50]

	return referenceMean


# def getContigsFromFile(contigFileName):
    # nameList = []
    # origSize = []
    # contigNum = 0
    # currentChrom = ""
    # seq = ""
    # contigDict = {}

    # try:
        # incontigfile = open(contigFileName)
    # except IOError:
        # print "Error opening contig file: %s" % contigFileName
        # return contigNum, nameList, contigDict, origSize

    # for line in incontigfile:
        # if ">" in line:
            # if currentChrom !="":
                # nameList.append(currentChrom)
                # contigDict[contigNum] = seq
                # origSize.append(len(seq))
                # contigNum += 1

            # currentChrom = line.strip().split()[0][1:]
            # seq = ""
        # else:
            # seq += line.strip()

    # incontigfile.close()

    # return contigNum, nameList, contigDict, origSize

# Here is the modified version of getContigsFromFile() added in RNAPATH*, original version fails to load the last sequence
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

	with open(contigFileName, 'r') as incontigfile:
		for line in incontigfile:
			if ">" in line:
				if "\r" in line:
					line = line.strip("\n")
				chrom = line[1:-1]
				nameList.append(chrom)
				contigNum += 1
				if seq :
					prevChrom = nameList[contigNum-2]
					contigDict[contigNum-2]=seq
					origSize.append(len(seq))
					seq=""
			else:
				seq += line.strip()

	contigDict[contigNum-1]=seq
	origSize.append(len(seq))
	incontigfile.close()
	return contigNum, nameList, contigDict, origSize

def getPath(contigNum, distalPairsfile, nameList):
	edgeMatrix = EdgeMatrix(contigNum)

	print len(edgeMatrix.edgeArray)
	try:
		print len(edgeMatrix.edgeArray[50])
	except IndexError:
		pass

	print "processing distal pairs"
	verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = processDistalPairsFile(distalPairsfile, edgeMatrix, nameList)

#print "notSoloDict",notSoloDict #added
#print "verticesWithEdges",verticesWithEdges #added

	willVisitList = verticesWithEdges.keys()
	willVisitList.sort() #willVisitList before [0, 1, 2, 3, 4, 5, 6]
	print "visiting %d vertices" % len(willVisitList)

	print "cleaning up graph of edges with weight 1"
	verticesToDelete = []
#print "edgeMatrix before", edgeMatrix.edgeArray #added

#print "vertexEdges",vertexEdges #added, vertexEdges {0: [3, 4, 5, 6], 1: [2], 2: [1], 3: [0], 4: [0], 5: [0], 6: [0]}

	for rindex in willVisitList: # if a contig is not in notSoloDict, it means that all connections to/from this contig have weight < 1(cutoff), and those weight will be converted to Zero in edgeMatrix.edgeArray
		if rindex not in notSoloDict:
			cindex = vertexEdges[rindex][0]
			edgeMatrix.edgeArray[rindex][cindex] = 0
			edgeMatrix.edgeArray[cindex][rindex] = 0
			verticesToDelete.append(rindex)
#print "edgeMatrix.edgeArray after 1", edgeMatrix.edgeArray


	for vertex in verticesToDelete:
		willVisitList.remove(vertex)   
		
	print "%d 1-edges zeroed out" % len(verticesToDelete)

	zeroedEdge = 0
	print "visiting %d vertices" % len(willVisitList)

	leafList = []
#print "willVisitList after",willVisitList
	print "picking top 2 edges per vertex - zero out others"
	for rindex in willVisitList:
		vertices = vertexEdges[rindex]
		rEdges = []
		for avertex in vertices:
			if avertex in willVisitList:
				rEdges.append((edgeMatrix.edgeArray[rindex][avertex], avertex))

		if len(rEdges) > 2:
			rEdges.sort()
			rEdges.reverse()
			zeroedEdge += len(rEdges[2:])
			for (weight, cindex) in rEdges[2:]:
				edgeMatrix.edgeArray[rindex][cindex] = 0  #further zero out the non-top2-weight edges
				edgeMatrix.edgeArray[cindex][rindex] = 0
			leafList.append(rindex) #added in RNAPATH*
		elif len(rEdges) == 1:
			if edgeMatrix.edgeArray[rindex][rEdges[0][1]] > 1:
				leafList.append(rindex)
		else:
			leafList.append(rindex) #added in RNAPATH*

			
#print "leafList", leafList #added
#print "edgeMatrix.edgeArray",edgeMatrix.edgeArray #added

	print "zeroed out %d lower-weight edges at vertices with degree > 2" % zeroedEdge
	pathList, visitedDict = traverseGraph(leafList, edgeMatrix)

	return pathList, edgeSenseDict, visitedDict

def traverseGraph(leafList, edgeMatrix):
	pathList = []
	visitedDict = {}
	leafList.sort()
	print "traveling through the graph"
	for rindex in leafList:
		if visitedDict.has_key(rindex):
			pass
		else:
			#path = edgeMatrix.visitLink(rindex)# orig
			path = edgeMatrix.visitLink(rindex,visitedDict.keys()) #added
			if len(path) > 1:
				for vertex in path:
					visitedDict[vertex] = ""

				print path
				pathList.append(path)

	return pathList, visitedDict

def processDistalPairsFile(distalPairsfilename, edgeMatrix, nameList):
	contigToRowLookup = {}
	verticesWithEdges = {}
	vertexEdges = {}
	notSoloDict = {}
	edgeSenseDict = {}

	distalPairs = open(distalPairsfilename)
	for line in distalPairs:
		if line[0] == "#":
			continue

		fields = line.strip().split()
		contA = "%s" % fields[1]
		try:
			contig1 = contigToRowLookup[contA]
		except KeyError:
			try:
				contig1 = nameList.index(contA)
				contigToRowLookup[contA] = contig1
			except ValueError:
				print "problem with end1: ", line
				continue

		sense1 = fields[3]

		contB = "%s" % fields[4]
		try:
			contig2 = contigToRowLookup[contB]
		except KeyError:
			try:
				contig2 = nameList.index(contB)
				contigToRowLookup[contB] = contig2
			except ValueError:
				print "problem with end2: ", line
				continue

		sense2 = fields[6]

		edgeMatrix.edgeArray[contig1][contig2] += 1
		edgeMatrix.edgeArray[contig2][contig1] += 1
		verticesWithEdges[contig1] = ""
		verticesWithEdges[contig2] = ""
		if (contig1, contig2) in edgeSenseDict:
			edgeSenseDict[contig1, contig2].append((sense1, sense2))
		elif (contig2, contig1) in edgeSenseDict:
			edgeSenseDict[contig2, contig1].append((sense2, sense1))
		else:
			edgeSenseDict[contig1, contig2] = [(sense1, sense2)]

		if contig1 in vertexEdges:
			if contig2 not in vertexEdges[contig1]:
				vertexEdges[contig1].append(contig2)
		else:
			vertexEdges[contig1] = [contig2]

		if contig2 in vertexEdges:
			if contig1 not in vertexEdges[contig2]:
				vertexEdges[contig2].append(contig1)
		else:
			vertexEdges[contig2] = [contig1]

		if edgeMatrix.edgeArray[contig1][contig2] > 1:
			notSoloDict[contig1] = "" #if a contig has at least one connections with weight >1, it will be added to the notSoloDict
									  # in other words, if all connections to/from a contig have weight <1, it will not be visited in scaffolding path
			notSoloDict[contig2] = ""

	distalPairs.close()

	return verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict

class EdgeMatrix:
	""" Describes a sparse matrix to hold edge data.
	"""

	def __init__(self, dimension):
		self.dimension = dimension
		self.edgeArray = zeros((self.dimension, self.dimension), int16)

	def visitLink(self, fromVertex, ignoreList=[]):
		returnPath = [fromVertex]
		toVertex = []
		for toindex in xrange(self.dimension):
			if self.edgeArray[fromVertex][toindex] > 1 and toindex not in ignoreList:
				toVertex.append(toindex)

		for vertex in toVertex:
			if sum(self.edgeArray[vertex]) == self.edgeArray[fromVertex][vertex]:
				self.edgeArray[fromVertex][vertex] = 0
				self.edgeArray[vertex][fromVertex] = 0
				return returnPath + [vertex]
			else:
				self.edgeArray[fromVertex][vertex] = 0
				try:
					return returnPath + self.visitLink(vertex, returnPath)
				except IOError:
					return returnPath + [vertex]
		return []

if __name__ == "__main__":
    main(sys.argv)
