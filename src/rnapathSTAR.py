import os
import sys
import math
import string
import collections
from numpy import zeros, int16

from lib import agouti_log as agLOG
from lib import agouti_gff as agGFF

class Graph:
	def __init__(self, dimension):
		self.dimension = dimension
		self.edgeArray = zeros((self.dimension, self.dimension), int16)

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

def scaffolding(walkStarts, graph):
	scafPaths = []
	visitedDict = {}
	walkStarts.sort()
	ignoreList = []
	for rindex in walkStarts:
		if visitedDict.has_key(rindex):
			pass
		else:
			moduleDEBUGLogger.debug(">graph walk start from node: %s" %(rindex))
			path, visitedDict = graph.walk_graph(rindex, visitedDict) #added
			if len(path) > 1:
				moduleDEBUGLogger.debug("return path: %s" %(",".join(map(str, path))))
				moduleDEBUGLogger.debug("pathLen=%d" %(len(path)))
				scafPaths.append(path)
	moduleDEBUGLogger.debug("number of paths: %d" %(len(scafPaths)))
	for k in visitedDict.iterkeys():
		moduleDEBUGLogger.debug("%s" %(k))
	moduleProgressLogger.info("number of visited nodes: %d" %(len(visitedDict)))
	return scafPaths, visitedDict

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

def reduce_complexity(willVisitList, vertexEdges, graph,
					  moduleOutDir, prefix,
					  seqNames, minSupport):
	moduleProgressLogger.info("Simplifying graph")
	moduleDEBUGLogger.debug("Simplifying graph")
	outGraphFile = os.path.join(moduleOutDir, "%s.rnapathSTAR.graph.gv" %(prefix))
	dGRAPH = {}
	zeroedEdge = 0
	walkStarts = []
	with open(outGraphFile, 'w') as fGRAPH:
		fGRAPH.write("graph {\n")
		for vertexA in willVisitList:
			vertices = vertexEdges[vertexA]
			rEdges = []
			for avertex in vertices:
				if avertex in willVisitList:
					rEdges.append((graph.edgeArray[vertexA][avertex], avertex))
			moduleDEBUGLogger.debug("vertexA %s - Edges %s" %(seqNames[vertexA], rEdges))

			verticesToDelete = []
			for i in range(len(rEdges)):
				weight, vertexB = rEdges[i]
				moduleDEBUGLogger.debug("\tvertexB %s - weight %d" %(seqNames[vertexB], weight))
				if i >= 3:
					graph.edgeArray[vertexA][vertexB] = 0
					graph.edgeArray[vertexB][vertexA] = 0
					moduleDEBUGLogger.debug("\t[reduce complexity] - remove edge %s - %s"
											%(seqNames[vertexA], seqNames[vertexB]))
					if (vertexA, vertexB) not in dGRAPH and \
					   (vertexB, vertexA) not in dGRAPH:
						dGRAPH[vertexA, vertexB] = 1
						fGRAPH.write("\t%s -- %s[label=%s, color=purple, penwidth=1];\n"
									 %(seqNames[vertexA], seqNames[vertexB],
									 graph.edgeArray[vertexA][vertexB]))
				else:
					if weight < minSupport:
						if (vertexA, vertexB) not in dGRAPH and \
						   (vertexB, vertexA) not in dGRAPH:
							dGRAPH[vertexA, vertexB] = 1
							fGRAPH.write("\t%s -- %s[label=%s, color=red, penwidth=1];\n"
										 %(seqNames[vertexA], seqNames[vertexB],
										 graph.edgeArray[vertexA][vertexB]))
						verticesToDelete.append(vertexB)
						graph.edgeArray[vertexA][vertexB] = 0
						graph.edgeArray[vertexB][vertexA] = 0
						moduleDEBUGLogger.debug("\t[insufficent support] - remove edge %s - %s"
												%(seqNames[vertexA], seqNames[vertexB]))
					else:
						if (vertexA, vertexB) not in dGRAPH and \
						   (vertexB, vertexA) not in dGRAPH:
							dGRAPH[vertexA, vertexB] = 1
							fGRAPH.write("\t%s -- %s[label=%s, color=blue, penwidth=2];\n"
										 %(seqNames[vertexA], seqNames[vertexB],
										 graph.edgeArray[vertexA][vertexB]))
			if len(rEdges)-len(verticesToDelete) > 0:
				walkStarts.append(vertexA)							#added in RNAPATH*
			if len(verticesToDelete) >= 1:
				for vertex in verticesToDelete:
					vertexEdges[vertexA].remove(vertex)
					vertexEdges[vertex].remove(vertexA)
		fGRAPH.write("}")
	moduleProgressLogger.info("Number of starting nodes: %d" %(len(walkStarts)))

	return walkStarts, graph

def build_graph(joinPairsFile, graph, seqNames):
	moduleProgressLogger.info("Building graph from joining reads pairs")

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

				graph.edgeArray[contig1][contig2] += 1
				graph.edgeArray[contig2][contig1] += 1
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

				if graph.edgeArray[contig1][contig2] > 1:
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
	moduleProgressLogger.info("Initializing edge-weighted graph")
	graph = Graph(nContig)
	moduleDEBUGLogger.debug("Dimension of edge matrix: %d x %d" %(nContig, nContig))

	verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = build_graph(joinPairsFile, graph, seqNames)

	willVisitList = verticesWithEdges.keys()
	willVisitList.sort()
	moduleProgressLogger.debug("%d vertices in the graph" %(len(willVisitList)))
	moduleDEBUGLogger.debug("%d vertices in the graph" %(len(willVisitList)))

	walkStarts, graph = reduce_complexity(willVisitList, vertexEdges, graph,
											   moduleOutDir, prefix,
											   seqNames, minSupport)

	moduleProgressLogger.info("Start graph walk")
	moduleDEBUGLogger.debug("Start graph walk")
	scafPaths, visitedDict = scaffolding(walkStarts, graph)

	moduleProgressLogger.info("%d paths scaffolded" %(len(scafPaths)))
	moduleProgressLogger.info("Succeeded")

	return scafPaths, edgeSenseDict, visitedDict
