import os
import sys
import math
import string
import collections
import itertools
import operator
#from numpy import zeros, int16

from lib import agouti_log as agLOG
from lib import agouti_gff as agGFF

class Graph(object):
	def __init__(self, outGraphFile, graph={}):
		self.outGraphFile = outGraphFile
		self.graph = graph
		self.weights = {}
		self.senses = {}

	def build_graph(self, joinPairsFile, vertex2Name):
		moduleProgressLogger.info("Building graph from joining reads pairs")
		with open(joinPairsFile, 'r') as fJOINPAIRS:
			for line in fJOINPAIRS:
				if not line.startswith('#'):
					tmp_line = line.strip().split()

					contigA = tmp_line[1]
					vertexA = vertex2Name.index(contigA)
					senseA = tmp_line[3]
					contigB = tmp_line[4]
					vertexB = vertex2Name.index(contigB)
					senseB = tmp_line[6]

					moduleDEBUGLogger.debug("vertexA %s | vertexB %s" %(contigA, contigB))
					self.add_vertices(vertexA, vertexB)
					self.add_edge(vertexA, vertexB)
					self.update_weight(vertexA, vertexB)
					self.add_sense(vertexA, senseA, vertexB, senseB)
					#moduleDEBUGLogger.debug("EdgeA %s" %(",".join([vertex2Name[k] for k in self.graph[vertexA]])))
					#moduleDEBUGLogger.debug("EdgeB %s" %(",".join([vertex2Name[k] for k in self.graph[vertexB]])))
					#moduleDEBUGLogger.debug("Weight %d" %(self.weights[vertexA, vertexB]))
					#moduleDEBUGLogger.debug("sense %s" %(self.senses[vertexA, vertexB]))

	def add_vertices(self, *vertices):
		for vertex in vertices:
			if vertex not in self.graph:
				self.graph[vertex] = []

	def add_edge(self, vertexA, vertexB):
		if (vertexA, vertexB) not in self.weights:
			self.graph[vertexA].append(vertexB)
			self.weights[vertexA, vertexB] = 0
		if (vertexB, vertexA) not in self.weights:
			self.graph[vertexB].append(vertexA)
			self.weights[vertexB, vertexA] = 0

	def update_weight(self, vertexA, vertexB):
		self.weights[vertexA, vertexB] += 1
		self.weights[vertexB, vertexA] += 1

	def add_sense(self, vertexA, senseA, vertexB, senseB):
		if (vertexA, vertexB) in self.senses:
			self.senses[vertexA, vertexB].append((senseA, senseB))
		elif (vertexB, vertexA) in self.senses:
			self.senses[vertexB, vertexA].append((senseB, senseA))
		else:
			self.senses[vertexA, vertexB] = [(senseA, senseB)]

	def remove_edge(self, vertexA, vertexB):
		self.graph[vertexA].remove(vertexB)
		self.graph[vertexB].remove(vertexA)
		del self.weights[vertexA, vertexB]
		del self.weights[vertexB, vertexA]
		if (vertexA, vertexB) in self.senses:
			del self.senses[vertexA, vertexB]
		elif (vertexB, vertexA) in self.senses:
			del self.senses[vertexB, vertexA]

	def get_vertices(self):
		return list(self.graph.iterkeys())

	def simplify(self, vertex2Name, minSupport):
		moduleProgressLogger.info("Simplifying graph")
		moduleDEBUGLogger.debug("Simplifying graph")
		self.leaves = []
		startVertices = []
		self.nonLeaves = []
		for vertexA in self.get_vertices():
			moduleDEBUGLogger.debug("vertexA - %s - %s" %(vertexA, vertex2Name[vertexA]))
			vertices = self.graph[vertexA]
			edges = []
			for vertexB in vertices:
				if vertexB in vertices:
					weight = self.weights[vertexA, vertexB]
					moduleDEBUGLogger.debug("----vertexB %s - %s - weight %d" %(vertexB, vertex2Name[vertexB], weight))
					edges.append((weight, vertexB))
			moduleDEBUGLogger.debug("----edges %s" %(str(edges)))
			if len(edges) >= 2:
				edges.sort(reverse=True)
				for (weight, vertexB) in edges[2:]:
					self.remove_edge(vertexA, vertexB)
				moduleDEBUGLogger.debug("----add to nonLeaves")
				self.nonLeaves.append(vertexA)
			elif len(edges) == 1:
				if self.weights[vertexA, edges[0][1]] >= minSupport:
					self.leaves.append(vertexA)
					startVertices.append(vertexA)
					moduleDEBUGLogger.debug("----add to leafs")

	def dfs(self, vertex, vertex2Name, minSupport, seen=None, path=None):
		if seen is None: seen = [vertex]
		if path is None: path = [vertex]

		for vertexB in self.graph[vertex]:
			weight = self.weights[vertex, vertexB]
			if weight >= minSupport and vertexB not in seen:
				seen.append(vertexB)
				path.append(vertexB)
				path = self.dfs(vertexB, vertex2Name, minSupport, seen, path)
		return path

	def graph2dot(self, scafPaths, vertex2Name, minSupport):
		with open(self.outGraphFile, 'w') as fGRAPH:
			dPairs = {}
			fGRAPH.write("graph {\n")
			for path in scafPaths:
				#fGRAPH.write("\t%s[color=red, penwidth=3.0]\n" %(" -- ".join([vertex2Name[k] for k in path])))
				for i in range(1, len(path)):
					vertexA = path[i-1]
					vertexB = path[i]
					weight = self.weights[vertexA, vertexB]
					dPairs[vertexA, vertexB] = 1
					dPairs[vertexB, vertexA] = 1
					fGRAPH.write("\t%s -- %s[label=\"%d\", weight=\"%d\", color=red, penwidth=3.0]\n"
								 %(vertex2Name[vertexA], vertex2Name[vertexB], weight, weight))

			for (vertexA, vertexB) in self.weights.iterkeys():
				weight = self.weights[vertexA, vertexB]
				if (vertexA, vertexB) not in dPairs:
					dPairs[vertexA, vertexB] = 1
					dPairs[vertexB, vertexA] = 1
					if weight < minSupport:
						fGRAPH.write("\t%s -- %s[style=\"dotted\", color=black, penwidth=1.0]\n"
									 %(vertex2Name[vertexA], vertex2Name[vertexB]))
					else:
						fGRAPH.write("\t%s -- %s[label=\"%d\", weight=\"%d\", color=black, penwidth=2.0]\n"
									 %(vertex2Name[vertexA], vertex2Name[vertexB], weight, weight))
			fGRAPH.write("}")

class RNAPATHSTAR_Graph(Graph):
	def start(self, joinPairsFile, vertex2Name, dCtgPair2GenePair, minSupport):
		self.build_graph(joinPairsFile, vertex2Name)
		vertices = self.get_vertices()
		moduleProgressLogger.info("%d vertices in the graph" %(len(vertices)))
		moduleDEBUGLogger.debug("%d vertices in the graph" %(len(vertices)))

		self.simplify(vertex2Name, minSupport)

		moduleProgressLogger.info("Start graph walk")
		scafPaths = self.scaffolding(vertex2Name, dCtgPair2GenePair, minSupport)
		return scafPaths

	def walk_graph(self, fromVertex, vertex2Name, minSupport, visitedDict={}):
		returnPath = [fromVertex]
		print "fromVertex", vertex2Name[fromVertex]
		print "returnPath", [vertex2Name[k] for k in returnPath]
		visitedDict[fromVertex] = ""
		toVertex = []
		for toindex in self.graph[fromVertex]:
			weight = self.weights[fromVertex, toindex]
			if weight >= minSupport and toindex not in visitedDict:
				toVertex.append((weight, toindex))

		toVertex.sort(reverse=True)
		print "toVertex", toVertex
		for (weight, vertex) in toVertex:
			if vertex in vertex2Name:
				print "vertex", vertex2Name[vertex]
			else:
				print "vertex", vertex
			totalWeight = 0
			for k in self.graph[vertex]:
				totalWeight += self.weights[vertex, k]
			if totalWeight == self.weights[fromVertex,vertex]:
				self.weights[fromVertex, vertex] = 0
				self.weights[vertex, fromVertex] = 0
				returnPath += [vertex]
				visitedDict[vertex] = ""
				return returnPath, visitedDict
			else:
				self.weights[fromVertex, vertex] = 0
				try:
					path, visitedDict = self.walk_graph(vertex, vertex2Name, minSupport, dict(visitedDict, **{vertex:""}))
					returnPath += path
					visitedDict[vertex] = ""
					return returnPath, visitedDict
				except IOError:
					returnPath += [vertex]
					return returnPath, visitedDict
		return returnPath, visitedDict

	def scaffolding(self, vertex2Name, dCtgPair2GenePair, minSupport):
		scafPaths = []
		dVisited = {}
		startVertices = self.leaves
		for vertex in startVertices:
			if dVisited.has_key(vertex):
				pass
			else:
				moduleDEBUGLogger.debug(">graph walk start from node: %s" %(vertex))
				path, dVisited = self.walk_graph(vertex, vertex2Name, minSupport, dVisited) #added
				if len(path) > 1:
					moduleDEBUGLogger.debug("return path: %s" %(",".join(map(str, path))))
					moduleDEBUGLogger.debug("pathLen=%d" %(len(path)))
					scafPaths.append(path)
		moduleDEBUGLogger.debug("number of paths: %d" %(len(scafPaths)))

		dSense = {}
		for (fVertex, tVertex), senses in self.senses.iteritems():
			counter = collections.Counter(senses)
			sense = sorted(counter.items(), key=operator.itemgetter(1), reverse=True)[0][0]
			dSense[fVertex, tVertex] = sense

		moduleDEBUGLogger.debug("Finding paths from non-leaves")
		leftOvers = [k for k in startVertices if k not in dVisited]
		leftOvers += [k for k in self.nonLeaves if k not in dVisited and k not in leftOvers]
#		leftOvers = [k for k in nonLeaves.iterkeys() if k not in dVisited]
#		leftOvers += [k for k in willVisitList if k not in dVisited and k not in leftOvers]
		moduleDEBUGLogger.debug("%d leftovers" %(len(leftOvers)))

		dBestPaths = {}
		visitedLeftOvers = []
		for vertex in leftOvers:
			if vertex not in dVisited:
				moduleDEBUGLogger.debug("%d - %s has not visited yet" %(vertex, vertex2Name[vertex]))
				path = self.dfs(vertex, vertex2Name, minSupport)
				moduleDEBUGLogger.debug("missing path: %s" %(",".join([vertex2Name[k] for k in path])))
				vertices = ""
				if len(path) == 1:
					dVisited[vertex] = ""
					visitedLeftOvers.append(vertex)
					moduleDEBUGLogger.debug("number of visited nodes: %d" %(len(dVisited)))
				elif len(path) > 1:
					moduleDEBUGLogger.debug("----before %d" %(len(dVisited)))
					for k in path:
						dVisited[k] = ""
					moduleDEBUGLogger.debug("----after %d" %(len(dVisited)))
					maxPathLen = 0
					bestPath = []
					for i in range(len(path)):
						if i > 0:
							tmpPath = self.dfs(path[i], vertex2Name, minSupport)
						else:
							tmpPath = path
						moduleDEBUGLogger.debug("walking path: %s" %(",".join([vertex2Name[k] for k in tmpPath])))

						curVertex = tmpPath[0]
						curCtg = vertex2Name[curVertex]
						curSensePair = ""
						possiblePath = []
						preSensePair = ""
						for i in range(1, len(tmpPath)):
							nextVertex = tmpPath[i]
							nextCtg = vertex2Name[nextVertex]
							moduleDEBUGLogger.debug("----cur %s | next %s" %(curCtg, nextCtg))
							if (curVertex, nextVertex) in dSense:
								ctgPair = (curVertex, nextVertex)
								genePair = dCtgPair2GenePair[ctgPair]
								curSensePair = "".join(dSense[ctgPair])
							elif (nextVertex, curVertex) in dSense:
								ctgPair = (nextVertex, curVertex)
								genePair = dCtgPair2GenePair[ctgPair][::-1]
								curSensePair = "".join(dSense[ctgPair][::-1])
								moduleDEBUGLogger.debug("----reverse contig pair")
							else:
								moduleDEBUGLogger.debug("----the pair has been removed, stop extension")
								break
							moduleDEBUGLogger.debug("----gene pair: %s %s" %(genePair[0].geneID, genePair[1].geneID))
							if preSensePair == "":
								preSensePair = curSensePair
								moduleDEBUGLogger.debug("----previous sense pair: %s" %(str(preSensePair)))
								moduleDEBUGLogger.debug("----current sense pair: %s" %(str(curSensePair)))
							else:
								moduleDEBUGLogger.debug("----previous sense pair: %s" %(str(preSensePair)))
								moduleDEBUGLogger.debug("----current sense pair: %s" %(str(curSensePair)))
								if preSensePair == "+-" and (curSensePair == "+-" or curSensePair == "++"):
									pass
								elif preSensePair == "++" and curSensePair == "--":
									pass
								elif preSensePair == "--" and curSensePair == "+-":
									pass
								elif preSensePair == "-+" and (curSensePair == "-+" or curSensePair == "--"):
									pass
								else:
									break
							possiblePath.append(curVertex)
							curVertex = nextVertex
							curCtg = vertex2Name[curVertex]
						possiblePath.append(curVertex)
						moduleDEBUGLogger.debug("----possiblePath: %s" %(",".join([vertex2Name[k] for k in possiblePath])))

						if len(possiblePath) == len(path):
							bestPath = possiblePath
							break
						else:
							if maxPathLen == 0:
								maxPathLen = len(possiblePath)
								bestPath = possiblePath
							else:
								if maxPathLen < len(possiblePath):
									maxPathLen = len(possiblePath)
									bestPath = possiblePath
						moduleDEBUGLogger.debug("----tmp bestPath: %s" %(",".join([vertex2Name[k] for k in bestPath])))
					moduleDEBUGLogger.debug("----Best path: %s" %(",".join([vertex2Name[k] for k in bestPath])))
					scafPaths.append(bestPath)
					moduleDEBUGLogger.debug("----before %d" %(len(dVisited)))
					for vertex in bestPath:
						visitedLeftOvers.append(k)
					moduleDEBUGLogger.debug("----before %d" %(len(dVisited)))
					moduleDEBUGLogger.debug("----after %d" %(len(dVisited)))

		moduleDEBUGLogger.debug("%d visited leftovers %s" %(len(visitedLeftOvers), str([vertex2Name[k] for k in visitedLeftOvers])))

		moduleDEBUGLogger.debug("number of visited nodes: %d" %(len(dVisited)))
		moduleProgressLogger.info("number of visited nodes: %d" %(len(dVisited)))
		moduleDEBUGLogger.debug("number of paths scaffolded: %d" %(len(scafPaths)))
		return scafPaths

class AGOUTI_Graph(Graph):
	def start(self, joinPairsFile, vertex2Name, dCtgPair2GenePair, minSupport):
		self.build_graph(joinPairsFile, vertex2Name)
		vertices = self.get_vertices()
		moduleProgressLogger.info("%d vertices in the graph" %(len(vertices)))
		moduleDEBUGLogger.debug("%d vertices in the graph" %(len(vertices)))

		self.simplify(vertex2Name, minSupport)
		startVertices = self.get_start_vertices()
		moduleProgressLogger.info("%d vertices to start" %(len(self.get_start_vertices())))
		moduleDEBUGLogger.debug("%d vertices to start" %(len(self.get_start_vertices())))

		moduleProgressLogger.info("Start graph walk")
		scafPaths = self.scaffolding(vertex2Name, dCtgPair2GenePair, minSupport)

		moduleProgressLogger.info("Generating graph in DOT format")
		if self.outGraphFile != "" and self.outGraphFile is not None:
			self.graph2dot(scafPaths, vertex2Name, minSupport)
		return scafPaths

	def get_start_vertices(self):
		self.startVertices = [k for k in self.leaves] + \
							 [k for k in self.nonLeaves if k not in self.leaves]
		return self.startVertices

	def scaffolding(self, vertex2Name, dCtgPair2GenePair, minSupport):
		scafPaths = []
		dVisited = {}

		dSense = {}
		for (vertexA, vertexB), senses in self.senses.iteritems():
			counter = collections.Counter(senses)
			sense = sorted(counter.items(), key=operator.itemgetter(1), reverse=True)[0][0]
			moduleDEBUGLogger.debug("vertexA %s | vertexB %s | %s" %(vertex2Name[vertexA], vertex2Name[vertexB], sense))
			dSense[vertexA, vertexB] = sense

		dBestPaths = {}
		for vertex in self.startVertices:
			if vertex not in dVisited:
				moduleDEBUGLogger.debug("%d - %s has not visited yet" %(vertex, vertex2Name[vertex]))
				path = self.dfs(vertex, vertex2Name, minSupport)
				moduleDEBUGLogger.debug("missing path: %s" %(",".join([vertex2Name[k] for k in path])))
				vertices = ""
				if len(path) == 1:
					dVisited[vertex] = ""
					moduleDEBUGLogger.debug("number of visited nodes: %d" %(len(dVisited)))
				elif len(path) > 1:
					moduleDEBUGLogger.debug("----before %d" %(len(dVisited)))
					for k in path:
						dVisited[k] = ""
					moduleDEBUGLogger.debug("----after %d" %(len(dVisited)))
					maxPathLen = 0
					bestPath = []
					for i in range(len(path)):
						if i > 0:
							tmpPath = self.dfs(path[i], vertex2Name, minSupport)
						else:
							tmpPath = path
						moduleDEBUGLogger.debug("walking path: %s" %(",".join([vertex2Name[k] for k in tmpPath])))

						curVertex = tmpPath[0]
						curCtg = vertex2Name[curVertex]
						curSensePair = ""
						possiblePath = []
						preSensePair = ""
						for i in range(1, len(tmpPath)):
							nextVertex = tmpPath[i]
							nextCtg = vertex2Name[nextVertex]
							moduleDEBUGLogger.debug("----cur %s | next %s" %(curCtg, nextCtg))
							if (curVertex, nextVertex) in dSense:
								ctgPair = (curVertex, nextVertex)
								genePair = dCtgPair2GenePair[ctgPair]
								curSensePair = "".join(dSense[ctgPair])
							elif (nextVertex, curVertex) in dSense:
								ctgPair = (nextVertex, curVertex)
								genePair = dCtgPair2GenePair[ctgPair][::-1]
								curSensePair = "".join(dSense[ctgPair][::-1])
								moduleDEBUGLogger.debug("----reverse contig pair")
							else:
								moduleDEBUGLogger.debug("----no gene pair between this contig pair, stop extension")
								break
							moduleDEBUGLogger.debug("----gene pair: %s %s" %(genePair[0].geneID, genePair[1].geneID))
							if preSensePair == "":
								preSensePair = curSensePair
								moduleDEBUGLogger.debug("----previous sense pair: %s" %(str(preSensePair)))
								moduleDEBUGLogger.debug("----current sense pair: %s" %(str(curSensePair)))
							else:
								moduleDEBUGLogger.debug("----previous sense pair: %s" %(str(preSensePair)))
								moduleDEBUGLogger.debug("----current sense pair: %s" %(str(curSensePair)))
								if preSensePair == "+-" and (curSensePair == "+-" or curSensePair == "++"):
									pass
								elif preSensePair == "++" and curSensePair == "--":
									pass
								elif preSensePair == "--" and curSensePair == "+-":
									pass
								elif preSensePair == "-+" and (curSensePair == "-+" or curSensePair == "--"):
									pass
								else:
									break
							possiblePath.append(curVertex)
							curVertex = nextVertex
							curCtg = vertex2Name[curVertex]
						possiblePath.append(curVertex)
						moduleDEBUGLogger.debug("----possiblePath: %s" %(",".join([vertex2Name[k] for k in possiblePath])))

						if len(possiblePath) == len(path):
							bestPath = possiblePath
							break
						else:
							if maxPathLen == 0:
								maxPathLen = len(possiblePath)
								bestPath = possiblePath
							else:
								if maxPathLen < len(possiblePath):
									maxPathLen = len(possiblePath)
									bestPath = possiblePath
						moduleDEBUGLogger.debug("----tmp bestPath: %s" %(",".join([vertex2Name[k] for k in bestPath])))
					moduleDEBUGLogger.debug("----Best path: %s" %(",".join([vertex2Name[k] for k in bestPath])))
					scafPaths.append(bestPath)

		moduleDEBUGLogger.debug("number of visited nodes: %d" %(len(dVisited)))
		moduleProgressLogger.info("number of visited nodes: %d" %(len(dVisited)))
		return scafPaths

def dfs(vertex, graph, vertex2Name, minSupport, seen=None, path=None):
	if seen is None: seen = [vertex]
	if path is None: path = [vertex]

	for t in xrange(graph.dimension):
		if graph.edgeArray[vertex][t] >= minSupport and t not in seen:
			seen.append(t)
			path.append(t)
			path = dfs(t, graph, vertex2Name, minSupport, seen, path)
	return path

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
					  vertex2Name, minSupport):
	moduleProgressLogger.info("Simplifying graph")
	moduleDEBUGLogger.debug("Simplifying graph")
	outGraphFile = os.path.join(moduleOutDir, "%s.rnapathSTAR.graph.gv" %(prefix))
	dGRAPH = {}
	zeroedEdge = 0
	startVertices = []
	nonLeaves = {}
	with open(outGraphFile, 'w') as fGRAPH:
		fGRAPH.write("graph {\n")
		for vertexA in willVisitList:
			moduleDEBUGLogger.debug("vertexA - %s - %s" %(vertexA, vertex2Name[vertexA]))
			vertices = vertexEdges[vertexA]
			rEdges = []
			for avertex in vertices:
				if avertex in willVisitList:
					weight = graph.edgeArray[vertexA][avertex]
					moduleDEBUGLogger.debug("----vertexB %s - %s - weight %d" %(avertex, vertex2Name[avertex], weight))
					rEdges.append((weight, avertex))

			if len(rEdges) >= 2:
				rEdges.sort(reverse=True)
				for (weight, vertexB) in rEdges[2:]:
					graph.edgeArray[vertexA][vertexB] = 0
					graph.edgeArray[vertexB][vertexA] = 0
				moduleDEBUGLogger.debug("----add to nonLeaves")
				nonLeaves[vertexA] = 1
			elif len(rEdges) == 1:
				if graph.edgeArray[vertexA][rEdges[0][1]] >= minSupport:
					startVertices.append(vertexA)
					moduleDEBUGLogger.debug("----add to leafs")
	moduleProgressLogger.info("Number of starting nodes: %d" %(len(startVertices)))

	return startVertices, graph, nonLeaves

def build_graph(joinPairsFile, graph, vertex2Name):
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

				contig1 = vertex2Name.index(tmp_line[1])
				sense1 = tmp_line[3]
				contig2 = vertex2Name.index(tmp_line[4])
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

def run_scaffolding(algorithm, vertex2Name, joinPairsFile,
				    dCtgPair2GenePair, moduleOutDir, prefix, minSupport):
	moduleProgressLogFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.progressMeter" %(prefix))
	moduleDebugLogFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.debug" %(prefix))
	moduleOutputFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.txt" %(prefix))
	global moduleProgressLogger
	moduleProgressLogger = agLOG.AGOUTI_LOG(moduleName).create_logger(moduleProgressLogFile)
	global moduleDEBUGLogger
	moduleDEBUGLogger = agLOG.AGOUTI_DEBUG_LOG(moduleName+"_DEBUG").create_logger(moduleDebugLogFile)
	moduleProgressLogger.info("[BEGIN] Scaffolding - Algorithm - %s priority" %(algorithm))

#	nContig = len(vertex2Name)
#	moduleProgressLogger.info("Initializing edge-weighted graph")
#	graph = Graph(nContig, minSupport)
#	moduleDEBUGLogger.debug("Dimension of edge matrix: %d x %d" %(nContig, nContig))

	#verticesWithEdges, vertexEdges, notSoloDict, edgeSenseDict = build_graph(joinPairsFile, graph, vertex2Name)
	outGraphFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.graph.dot" %(prefix))
	scafPaths = []
	if algorithm == "gene":
		graph = AGOUTI_Graph(outGraphFile)
		scafPaths = graph.start(joinPairsFile, vertex2Name, dCtgPair2GenePair, minSupport)
	elif algorithm == "weight":
		graph = RNAPATHSTAR_Graph(outGraphFile)
		scafPaths = graph.start(joinPairsFile, vertex2Name, dCtgPair2GenePair, minSupport)

#	willVisitList = verticesWithEdges.keys()
#	willVisitList.sort()
#	moduleProgressLogger.info("%d vertices in the graph" %(len(willVisitList)))
#	moduleDEBUGLogger.debug("%d vertices in the graph" %(len(willVisitList)))

#	startVertices, graph, nonLeaves = reduce_complexity(willVisitList, vertexEdges, graph,
#												  moduleOutDir, prefix,
#												  vertex2Name, minSupport)

#	moduleProgressLogger.info("Start graph walk")
#	moduleDEBUGLogger.debug("Start graph walk")
#	scafPaths, visitedDict = scaffolding(startVertices, graph, nonLeaves, vertex2Name, willVisitList, edgeSenseDict, dCtgPair2GenePair, minSupport)

	moduleProgressLogger.info("%d paths scaffolded" %(len(scafPaths)))
	moduleProgressLogger.info("Succeeded")

	return scafPaths, graph.senses
