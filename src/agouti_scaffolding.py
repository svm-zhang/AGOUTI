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
			if sum(self.edgeArray[vertex]) == self.edgeArray[fromVertex][vertex]:
				self.edgeArray[fromVertex][vertex] = 0
				self.edgeArray[vertex][fromVertex] = 0
				return returnPath + [vertex]
			else:
				self.edgeArray[fromVertex][vertex] = 0
				try:
					return returnPath + self.visitLink(vertex, ignoreList)
				except IOError:
					return returnPath + [vertex]
		return []

#	def sz_visitLink(self, fromVertex, ignoreList=[]):
	def sz_visitLink(self, fromVertex, visitedDict={}):
		returnPath = [fromVertex]
		print "fromVertex", fromVertex,
		print "returnPath", returnPath
		visitedDict[fromVertex] = ""
#		newVisited = [fromVertex]
		print "newVisited", visitedDict
		toVertex = []
		for toindex in xrange(self.dimension):
#			if self.edgeArray[fromVertex][toindex] > 1 and toindex not in ignoreList:
			if self.edgeArray[fromVertex][toindex] > 1 and toindex not in visitedDict:
				toVertex.append(toindex)
		print "toVertex", toVertex

		for vertex in toVertex:
			print "vertex", vertex
#			if vertex not in ignoreList:
			if sum(self.edgeArray[vertex]) == self.edgeArray[fromVertex][vertex]:
				self.edgeArray[fromVertex][vertex] = 0
				self.edgeArray[vertex][fromVertex] = 0
				returnPath += [vertex]
#				newVisited += [vertex]
				visitedDict[vertex] = ""
				return returnPath, visitedDict
#				return returnPath, newVisited
			else:
				self.edgeArray[fromVertex][vertex] = 0
				try:
					path, visited = self.sz_visitLink(vertex, dict(visitedDict, **{vertex:""}))
#					path, visited = self.sz_visitLink(vertex, newVisited+ignoreList)
					returnPath += path
#					newVisited += visited
					visitedDict[vertex] = ""
					return returnPath, visitedDict
#					return returnPath, newVisited
				except IOError:
					returnPath += [vertex]
					return returnPath, visitedDict
#					return returnPath, newVisited
		return returnPath, visitedDict
#		return returnPath, newVisited

def sz_traverseGraph(leafList, edgeMatrix):
	pathList = []
	visitedDict = {}
	leafList.sort()
	ignoreList = []
	print "traveling through the graph"
	for rindex in leafList:
		if visitedDict.has_key(rindex):
			pass
		else:
			print ">rindex", rindex
#			path = edgeMatrix.visitLink(rindex, visitedDict.keys()) # orig
			path, ignoreList = edgeMatrix.sz_visitLink(rindex, visitedDict.keys()) #added
			print "path", path, "ignore", ignoreList
			if len(path) > 1:
				visitedDict = dict(visitedDict, **{x: "" for x in ignoreList})
#				for vertex in path:
#					visitedDict[vertex] = ""
				print "path", path, "visited nodes", visitedDict.keys()
				pathList.append(path)
	print "number of visitedVertex: %d" %(len(visitedDict))
	print "number of paths: %d" %(len(pathList))
	sys.exit()
	return pathList, visitedDict

def traverseGraph(leafList, edgeMatrix):
	pathList = []
	visitedDict = {}
	leafList.sort()
	ignoreList = []
	print "traveling through the graph"
	for rindex in leafList:
		if visitedDict.has_key(rindex):
			pass
		else:
			print ">rindex", rindex
#			path = edgeMatrix.visitLink(rindex, visitedDict.keys()) # orig
			path, visitedDict = edgeMatrix.sz_visitLink(rindex, visitedDict) #added
#			path, ignoreList = edgeMatrix.sz_visitLink(rindex, visitedDict.keys()) #added
			if len(path) > 1:
#				visitedDict = dict(visitedDict, **{x: "" for x in ignoreList})
#				for vertex in path:
#					visitedDict[vertex] = ""
				print "path", path, "visited nodes", visitedDict.keys()
				pathList.append(path)
	print "number of visitedVertex: %d" %(len(visitedDict))
	print "number of paths: %d" %(len(pathList))
	return pathList, visitedDict

def check_orientation_conflicts(vertexA, vertexB, edgeSenseDict):
	fr, ff, rr, rf = 0, 0, 0, 0
	if (vertexA, vertexB) in edgeSenseDict:
		ff += edgeSenseDict[vertexA, vertexB].count(('+', '+'))
		fr += edgeSenseDict[vertexA, vertexB].count(('+', '-'))
		rr += edgeSenseDict[vertexA, vertexB].count(('-', '-'))
		rf += edgeSenseDict[vertexA, vertexB].count(('-', '+'))
	elif (vertexB, vertexA) in edgeSenseDict:
		ff += edgeSenseDict[vertexB, vertexA].count(('+', '+'))
		fr += edgeSenseDict[vertexB, vertexA].count(('+', '-'))
		rr += edgeSenseDict[vertexB, vertexA].count(('-', '-'))
		rf += edgeSenseDict[vertexB, vertexA].count(('-', '+'))
	counts = [fr, ff, rr, rf]
	print "vertexA", vertexA, "vertexB", vertexB, "counts", counts
	counts.sort(reverse=True)
	print "vertexA", vertexA, "vertexB", vertexB, "counts", counts
	if float(counts[1])/float(counts[0]) > 0.3:
		return True
	else:
		return False

def getPath(contigNum, joinPairsFile, nameList):
	edgeMatrix = EdgeMatrix(contigNum)

	print len(edgeMatrix.edgeArray)

	print "processing distal pairs"
	verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = process_join_pairs_file(joinPairsFile, edgeMatrix, nameList)

	for vertexA, vertices in vertexEdges.items():
		print "vertices", vertices
		toDelete = []
		for i in range(len(vertices)):
			vertexB = vertices[i]
			print "vertexA", vertexA, "vertexB", vertexB
			orientationConflict = check_orientation_conflicts(vertexA, vertexB, edgeSenseDict)
			if orientationConflict:
				toDelete.append(vertexB)
				edgeMatrix.edgeArray[vertexA][vertexB] = 0
				edgeMatrix.edgeArray[vertexB][vertexA] = 0
		for vertex in toDelete:
			vertexEdges[vertexA].remove(vertex)
			vertexEdges[vertex].remove(vertexA)

	willVisitList = verticesWithEdges.keys()
	willVisitList.sort()
	print "visiting %d vertices" % len(willVisitList)

	print "cleaning up graph of edges with weight 1"
	verticesToDelete = []

	for rindex in willVisitList: # if a contig is not in notSoloDict, it means that all connections to/from this contig have weight < 1(cutoff), and those weight will be converted to Zero in edgeMatrix.edgeArray
		if rindex not in notSoloDict:
			cindex = vertexEdges[rindex][0]
			edgeMatrix.edgeArray[rindex][cindex] = 0
			edgeMatrix.edgeArray[cindex][rindex] = 0
			verticesToDelete.append(rindex)
	print "solo to delete", len(verticesToDelete)

	for vertex in verticesToDelete:
		willVisitList.remove(vertex)
	print "%d 1-edges zeroed out" % len(verticesToDelete)

	zeroedEdge = 0
	print "visiting %d vertices" % len(willVisitList)

	leafList = []
	multi = []
	twos = []
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
				edgeMatrix.edgeArray[rindex][cindex] = 0	#further zero out the non-top2-weight edges
				edgeMatrix.edgeArray[cindex][rindex] = 0
			leafList.append(rindex)							#added in RNAPATH*
			multi.append(rindex)
		elif len(rEdges) == 1:
#			if edgeMatrix.edgeArray[rindex][rEdges[0][1]] > 1:
			leafList.append(rindex)
		else:
			twos.append(rindex)
			leafList.append(rindex) #added in RNAPATH*
	print "length leaf", len(leafList)
	print "length twos", len(twos)
	print "length multi", len(multi)
	print "twos", twos
	print "multi", multi
#		if rindex == 4222 or rindex == 4223:
#			print "vertices", vertices
#			print "orientationConflict", orientationConflict
#			print "rEdges", rEdges
#			print "leafList", leafList
#	sys.exit()

	print "zeroed out %d lower-weight edges at vertices with degree > 2" % zeroedEdge
	pathList, visitedDict = traverseGraph(leafList, edgeMatrix)

	return pathList, edgeSenseDict, visitedDict

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
#	for k, v in edgeSenseDict.items():
#		contig1 = k[0]
#		contig2 = k[1]

	return verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict

def agouti_scaffolding(contigNum, nameList, joinPairsFile):
	sys.stderr.write("Scaffolding ... \n")
	pathList, edgeSenseDict, visitedDict = getPath(contigNum, joinPairsFile, nameList)

	print pathList
	print "found %d paths" % len(pathList)

	return pathList, edgeSenseDict, visitedDict
