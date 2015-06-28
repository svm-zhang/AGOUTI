import os
import sys
import math
import string
import collections
from numpy import zeros, int16

from lib import agouti_gff as agGFF

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
			# this is the base case where no further extension can be found
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

def getPath(contigNum, joinPairsFile, nameList):
	edgeMatrix = EdgeMatrix(contigNum)

	print len(edgeMatrix.edgeArray)
#	try:
#		print len(edgeMatrix.edgeArray[50])
#	except IndexError:
#		pass

	print "processing distal pairs"
	verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = process_join_pairs_file(joinPairsFile, edgeMatrix, nameList)

	willVisitList = verticesWithEdges.keys()
	willVisitList.sort() #willVisitList before [0, 1, 2, 3, 4, 5, 6]
	print "visiting %d vertices" % len(willVisitList)

	print "cleaning up graph of edges with weight 1"
	verticesToDelete = []
#	print "edgeMatrix before", edgeMatrix.edgeArray #added

#	print "vertexEdges",vertexEdges #added, vertexEdges {0: [3, 4, 5, 6], 1: [2], 2: [1], 3: [0], 4: [0], 5: [0], 6: [0]}

	for rindex in willVisitList: # if a contig is not in notSoloDict, it means that all connections to/from this contig have weight < 1(cutoff), and those weight will be converted to Zero in edgeMatrix.edgeArray
		if rindex not in notSoloDict:
			cindex = vertexEdges[rindex][0]
			edgeMatrix.edgeArray[rindex][cindex] = 0
			edgeMatrix.edgeArray[cindex][rindex] = 0
			verticesToDelete.append(rindex)
			print nameList[rindex], nameList[cindex]
#	print "edgeMatrix.edgeArray after 1", edgeMatrix.edgeArray

	for vertex in verticesToDelete:
		willVisitList.remove(vertex)
	print "%d 1-edges zeroed out" % len(verticesToDelete)

	zeroedEdge = 0
	print "visiting %d vertices" % len(willVisitList)

	leafList = []
#	print "willVisitList after",willVisitList
	print "picking top 2 edges per vertex - zero out others"
	for rindex in willVisitList:
		vertices = vertexEdges[rindex]
		rEdges = []
		for avertex in vertices:
			if avertex in willVisitList:
				rEdges.append((edgeMatrix.edgeArray[rindex][avertex], avertex))

		if len(rEdges) > 2:
			rEdges.sort(reverse=True)
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
#			print nameList[rindex], nameList[rEdges[0][1]], nameList[rEdges[1][1]]
#	sys.exit()

#	print "leafList", leafList #added
#	print "edgeMatrix.edgeArray",edgeMatrix.edgeArray #added

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
#			path = edgeMatrix.visitLink(rindex) # orig
			path = edgeMatrix.visitLink(rindex,visitedDict.keys()) #added
			if len(path) > 1:
				for vertex in path:
					visitedDict[vertex] = ""
				print path
				pathList.append(path)
	return pathList, visitedDict

def process_join_pairs_file(joinPairsFilename, edgeMatrix, nameList):
	contigToRowLookup = {}
	verticesWithEdges = {}
	vertexEdges = {}
	notSoloDict = {}
	edgeSenseDict = {}

	fJOINPAIRS = open(joinPairsFilename)
	for line in fJOINPAIRS:
		if not line.startswith('#'):
			tmp_line = line.strip().split()

			contig1 = nameList.index(tmp_line[1])
			sense1 = tmp_line[3]
			contig2 = nameList.index(tmp_line[4])
			sense2 = tmp_line[6]

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

	fJOINPAIRS.close()
	for k, v in edgeSenseDict.items():
		contig1 = k[0]
		contig2 = k[1]
		if contig1 in [177, 5428] and \
		   contig2 in [177, 5428]:
			   print "177, 5428"
			   print k, v, "\n"
			   print "fuck"
		elif contig1 in [178, 5428] and \
			 contig2 in [178, 5428]:
				print "178, 5428"
				print k, v

	return verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict

def agouti_scaffolding(contigNum, nameList, joinPairsFile):
	sys.stderr.write("Scaffolding ... \n")
	pathList, edgeSenseDict, visitedDict = getPath(contigNum, joinPairsFile, nameList)

	print pathList
	print "found %d paths" % len(pathList)

	return pathList, edgeSenseDict, visitedDict
