import os
import sys
import math
import string
import collections
from numpy import zeros, int16

from lib import agouti_log as agLOG
from lib import agouti_gff as agGFF

class EdgeMatrix:
	"""
		Describes a sparse matrix to hold edge data.
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

	def walk_graph(self, fromVertex, visitedDict={}):
		returnPath = [fromVertex]
		print "fromVertex", fromVertex,
		print "returnPath", returnPath
		visitedDict[fromVertex] = ""
		toVertex = []
		for toindex in xrange(self.dimension):
			if self.edgeArray[fromVertex][toindex] > 1 and toindex not in visitedDict:
				toVertex.append(toindex)
		print "toVertex", toVertex

		for vertex in toVertex:
			print "vertex", vertex
			if sum(self.edgeArray[vertex]) == self.edgeArray[fromVertex][vertex]:
				self.edgeArray[fromVertex][vertex] = 0
				self.edgeArray[vertex][fromVertex] = 0
				returnPath += [vertex]
				visitedDict[vertex] = ""
				return returnPath, visitedDict
			else:
				self.edgeArray[fromVertex][vertex] = 0
				try:
					path, visitedDict = self.walk_graph(vertex, dict(visitedDict, **{vertex:""}))
					returnPath += path
					visitedDict[vertex] = ""
					return returnPath, visitedDict
				except IOError:
					returnPath += [vertex]
					return returnPath, visitedDict
		return returnPath, visitedDict

def sz_scaffolding(leafList, edgeMatrix):
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
			path, ignoreList = edgeMatrix.walk_graph(rindex, visitedDict.keys()) #added
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

def scaffolding(leafList, edgeMatrix):
	pathList = []
	visitedDict = {}
	leafList.sort()
	ignoreList = []
	for rindex in leafList:
		if visitedDict.has_key(rindex):
			pass
		else:
			moduleDEBUGLogger.debug(">graph walk start from node: %s" %(rindex))
			path, visitedDict = edgeMatrix.walk_graph(rindex, visitedDict) #added
			if len(path) > 1:
				moduleDEBUGLogger.debug("return path: %s" %(",".join(map(str, path))))
				pathList.append(path)
	moduleDEBUGLogger.debug("number of paths: %d" %(len(pathList)))
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

def getPath(nContig, joinPairsFile, seqNames, minSupport):
	moduleProgressLogger.info("Initializing edge-weighted graph")
	edgeMatrix = EdgeMatrix(nContig)

	moduleDEBUGLogger.debug("Dimension of edge matrix: %d x %d" %(nContig, nContig))

	verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = build_graph(joinPairsFile, edgeMatrix, seqNames)

	moduleProgressLogger.info("Simplifying graph")
	moduleDEBUGLogger.debug("Simplifying graph")
	for vertexA, vertices in vertexEdges.items():
		moduleDEBUGLogger.debug(">vertexA %s - connecting to %s" %(vertexA, ",".join(map(str, vertices))))
		toDelete = []
		for i in range(len(vertices)):
			vertexB = vertices[i]
			moduleDEBUGLogger.debug("\tvertexB %s - weight %d" %(vertexB, edgeMatrix.edgeArray[vertexA][vertexB]))
			if edgeMatrix.edgeArray[vertexA][vertexB] < minSupport:
				toDelete.append(vertexB)
				edgeMatrix.edgeArray[vertexA][vertexB] = 0
				edgeMatrix.edgeArray[vertexB][vertexA] = 0
				moduleDEBUGLogger.debug("\tDelete because of weight")
				continue
#			orientationConflict = check_orientation_conflicts(vertexA, vertexB, edgeSenseDict)
		if len(toDelete) > 0:
			moduleDEBUGLogger.debug("\tvertexA %s - delete connections %s" %(vertexA, ",".join(map(str, toDelete))))
		for vertex in toDelete:
			vertexEdges[vertexA].remove(vertex)
			vertexEdges[vertex].remove(vertexA)

	willVisitList = verticesWithEdges.keys()
	willVisitList.sort()
	moduleProgressLogger.debug("%d vertices in the graph" %(len(willVisitList)))
	moduleDEBUGLogger.debug("%d vertices in the graph" %(len(willVisitList)))

	verticesToDelete = []

	moduleProgressLogger.info("Preparing nodes to start graph walk")
	moduleDEBUGLogger.debug("Preparing nodes to start graph walk")
	zeroedEdge = 0
	leafList = []
	for rindex in willVisitList:
		vertices = vertexEdges[rindex]
		rEdges = []
		for avertex in vertices:
			if avertex in willVisitList:
				rEdges.append((edgeMatrix.edgeArray[rindex][avertex], avertex))

		moduleDEBUGLogger.debug("vertexA %s - Edges %s" %(rindex, rEdges))

		if len(rEdges) >= 2:
			rEdges.sort(reverse=True)
			zeroedEdge += len(rEdges[3:])
			for (weight, cindex) in rEdges[3:]:
				edgeMatrix.edgeArray[rindex][cindex] = 0
				edgeMatrix.edgeArray[cindex][rindex] = 0
			leafList.append(rindex)							#added in RNAPATH*
		elif len(rEdges) == 1:
			leafList.append(rindex)

	moduleProgressLogger.info("Start graph walk")
	moduleDEBUGLogger.debug("Start graph walk")
	pathList, visitedDict = scaffolding(leafList, edgeMatrix)

	return pathList, edgeSenseDict, visitedDict

def build_graph(joinPairsFile, edgeMatrix, seqNames):
	moduleProgressLogger.info("Building graph using joining reads pairs")

	contigToRowLookup = {}
	verticesWithEdges = {}
	vertexEdges = {}
	notSoloDict = {}
	edgeSenseDict = {}

	with open(joinPairsFile, 'r') as fJOINPAIRS:
		for line in fJOINPAIRS:
			if not line.startswith('#'):
				tmp_line = line.strip().split()

				contig1 = seqNames.index(tmp_line[1])
				sense1 = tmp_line[3]
				contig2 = seqNames.index(tmp_line[4])
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
					notSoloDict[contig1] = ""
					notSoloDict[contig2] = ""

	return verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict

def set_module_name(name):
	global moduleName
	moduleName = name

def rnapathSTAR(seqNames, joinPairsFile, moduleOutDir, prefix, minSupport):
	moduleProgressLogFile = os.path.join(moduleOutDir, "%s.rnapathSTAR.progressMeter" %(prefix))
	moduleDebugLogFile = os.path.join(moduleOutDir, "%s.rnapathSTAR.debug" %(prefix))
	moduleOutputFile = os.path.join(moduleOutDir, "%s.rnapathSTAR.txt" %(prefix))
	global moduleProgressLogger
	moduleProgressLogger = agLOG.AGOUTI_LOG(moduleName).create_logger(moduleProgressLogFile)
	global moduleDEBUGLogger
	moduleDEBUGLogger = agLOG.AGOUTI_DEBUG_LOG(moduleName+"_DEBUG").create_logger(moduleDebugLogFile)
	moduleProgressLogger.info("[BEGIN] Scaffolding")

	nContig = len(seqNames)
	pathList, edgeSenseDict, visitedDict = getPath(nContig, joinPairsFile, seqNames, minSupport)

#	moduleDEBUGLogger.debug(pathList)
	moduleProgressLogger.info("%d paths scaffolded" %(len(pathList)))
	moduleProgressLogger.info("Succeeded")

	return pathList, edgeSenseDict, visitedDict
