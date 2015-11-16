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

	def start_logger(self, moduleName, moduleOutDir, prefix, debug=0):
		"""
			initiate a logger to keep track stuff
		"""
		self.debug = debug
		progressLogFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.progressMeter" %(prefix))
		self.agSCAFProgress = agLOG.PROGRESS_METER(moduleName)
		self.agSCAFProgress.add_file_handler(progressLogFile)

		self.debugLogFile = None
		if self.debug:
			self.debugLogFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.debug" %(prefix))

	def build_graph(self, joinPairsFile, vertex2Name):
		"""
			build graph from given joining-pairs
		"""
		self.agSCAFProgress.logger.info("Building graph from joining reads pairs")
		if self.debug:
			buildDebug = agLOG.DEBUG("BUILD", self.debugLogFile)
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

					if (vertexA, vertexB) not in self.weights:
						buildDebug.debugger.debug("Building vertexA %s | vertexB %s" %(contigA, contigB))
					self.add_vertices(vertexA, vertexB)
					self.add_edge(vertexA, vertexB)
					self.update_weight(vertexA, vertexB)
					self.add_sense(vertexA, senseA, vertexB, senseB)
					#agSCAFDebug.debugger.debug("EdgeA %s" %(",".join([vertex2Name[k] for k in self.graph[vertexA]])))
					#agSCAFDebug.debugger.debug("EdgeB %s" %(",".join([vertex2Name[k] for k in self.graph[vertexB]])))
					#agSCAFDebug.debugger.debug("Weight %d" %(self.weights[vertexA, vertexB]))
					#agSCAFDebug.debugger.debug("sense %s" %(self.senses[vertexA, vertexB]))

	def add_vertices(self, *vertices):
		"""
			add the give number of contigs
			to the graph
		"""
		for vertex in vertices:
			if vertex not in self.graph:
				self.graph[vertex] = []

	def add_edge(self, vertexA, vertexB):
		"""
			add an edge between the give pair of contig
		"""
		if (vertexA, vertexB) not in self.weights:
			self.graph[vertexA].append(vertexB)
			self.weights[vertexA, vertexB] = 0
		if (vertexB, vertexA) not in self.weights:
			self.graph[vertexB].append(vertexA)
			self.weights[vertexB, vertexA] = 0

	def update_weight(self, vertexA, vertexB):
		"""
			add one weight to the give contig pair
		"""
		self.weights[vertexA, vertexB] += 1
		self.weights[vertexB, vertexA] += 1

	def add_sense(self, vertexA, senseA, vertexB, senseB):
		"""
			add sense info for the give pair of contig
		"""
		if (vertexA, vertexB) in self.senses:
			self.senses[vertexA, vertexB].append((senseA, senseB))
		elif (vertexB, vertexA) in self.senses:
			self.senses[vertexB, vertexA].append((senseB, senseA))
		else:
			self.senses[vertexA, vertexB] = [(senseA, senseB)]

	def remove_edge(self, vertexA, vertexB):
		"""
			remove an edge from the graph
		"""
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
		"""
			Simplify the graph by remove edges
			with low weights
		"""
		self.agSCAFProgress.logger.info("Simplifying graph")
		if self.debug:
			simplifyDebug = agLOG.DEBUG("SIMPLIFY", self.debugLogFile, 'a')
		self.leaves = []
		startVertices = []
		self.nonLeaves = []
		nEdgeRemoved = 0
		for vertexA in self.get_vertices():
			if self.debug:
				simplifyDebug.debugger.debug("vertexA - %s - neighbors - %s" %(vertexA, vertex2Name[vertexA]))
			vertices = self.graph[vertexA]
			edges = []
			for vertexB in vertices:
				if vertexB in vertices:
					weight = self.weights[vertexA, vertexB]
					if self.debug:
						simplifyDebug.debugger.debug("\tNeighbor B %s - %s - weight %d" %(vertexB, vertex2Name[vertexB], weight))
					edges.append((weight, vertexB))
			if len(edges) >= 2:
				edges.sort(reverse=True)
				nEdgeRemoved += len(edges[2:])
				for (weight, vertexB) in edges[2:]:
					self.remove_edge(vertexA, vertexB)
				if self.debug:
					simplifyDebug.debugger.debug("\tNon-leaf: Neighbor B %s - %s" %(vertexB, vertex2Name[vertexB]))
				self.nonLeaves.append(vertexA)
			elif len(edges) == 1:
				if self.weights[vertexA, edges[0][1]] >= minSupport:
					self.leaves.append(vertexA)
					startVertices.append(vertexA)
					if self.debug:
						simplifyDebug.debugger.debug("\tLeaf: Neighbor B %s - %s" %(vertexB, vertex2Name[vertexB]))
		self.agSCAFProgress.logger.info("%d Edges removed due to insufficient supports")

	def dfs(self, vertex, vertex2Name, minSupport, seen=None, path=None):
		"""
			Depth-first traversal from a given node
		"""
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
		"""
			output Graph in DOT format
		"""
		with open(self.outGraphFile, 'w') as fGRAPH:
			dPairs = {}
			fGRAPH.write("graph {\n")
			for path in scafPaths:
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
		agSCAFProgress.logger.info("%d vertices in the graph" %(len(vertices)))
		agSCAFDebug.debugger.debug("%d vertices in the graph" %(len(vertices)))

		self.simplify(vertex2Name, minSupport)

		agSCAFProgress.logger.info("Start graph walk")
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
				agSCAFDebug.debugger.debug(">graph walk start from node: %s" %(vertex))
				path, dVisited = self.walk_graph(vertex, vertex2Name, minSupport, dVisited) #added
				if len(path) > 1:
					agSCAFDebug.debugger.debug("return path: %s" %(",".join(map(str, path))))
					agSCAFDebug.debugger.debug("pathLen=%d" %(len(path)))
					scafPaths.append(path)
		agSCAFDebug.debugger.debug("number of paths: %d" %(len(scafPaths)))

		dSense = {}
		for (fVertex, tVertex), senses in self.senses.iteritems():
			counter = collections.Counter(senses)
			sense = sorted(counter.items(), key=operator.itemgetter(1), reverse=True)[0][0]
			dSense[fVertex, tVertex] = sense

		agSCAFDebug.debugger.debug("Finding paths from non-leaves")
		leftOvers = [k for k in startVertices if k not in dVisited]
		leftOvers += [k for k in self.nonLeaves if k not in dVisited and k not in leftOvers]
		agSCAFDebug.debugger.debug("%d leftovers" %(len(leftOvers)))

		dBestPaths = {}
		visitedLeftOvers = []
		for vertex in leftOvers:
			if vertex not in dVisited:
				agSCAFDebug.debugger.debug("%d - %s has not visited yet" %(vertex, vertex2Name[vertex]))
				path = self.dfs(vertex, vertex2Name, minSupport)
				agSCAFDebug.debugger.debug("missing path: %s" %(",".join([vertex2Name[k] for k in path])))
				vertices = ""
				if len(path) == 1:
					dVisited[vertex] = ""
					visitedLeftOvers.append(vertex)
					agSCAFDebug.debugger.debug("number of visited nodes: %d" %(len(dVisited)))
				elif len(path) > 1:
					agSCAFDebug.debugger.debug("\tbefore %d" %(len(dVisited)))
					for k in path:
						dVisited[k] = ""
					agSCAFDebug.debugger.debug("\tafter %d" %(len(dVisited)))
					maxPathLen = 0
					bestPath = []
					for i in range(len(path)):
						if i > 0:
							tmpPath = self.dfs(path[i], vertex2Name, minSupport)
						else:
							tmpPath = path
						agSCAFDebug.debugger.debug("walking path: %s" %(",".join([vertex2Name[k] for k in tmpPath])))

						curVertex = tmpPath[0]
						curCtg = vertex2Name[curVertex]
						curSensePair = ""
						possiblePath = []
						preSensePair = ""
						for i in range(1, len(tmpPath)):
							nextVertex = tmpPath[i]
							nextCtg = vertex2Name[nextVertex]
							agSCAFDebug.debugger.debug("\tcur %s | next %s" %(curCtg, nextCtg))
							if (curVertex, nextVertex) in dSense:
								ctgPair = (curVertex, nextVertex)
								genePair = dCtgPair2GenePair[ctgPair]
								curSensePair = "".join(dSense[ctgPair])
							elif (nextVertex, curVertex) in dSense:
								ctgPair = (nextVertex, curVertex)
								genePair = dCtgPair2GenePair[ctgPair][::-1]
								curSensePair = "".join(dSense[ctgPair][::-1])
								agSCAFDebug.debugger.debug("\treverse contig pair")
							else:
								agSCAFDebug.debugger.debug("\tthe pair has been removed, stop extension")
								break
							agSCAFDebug.debugger.debug("\tgene pair: %s %s" %(genePair[0].geneID, genePair[1].geneID))
							if preSensePair == "":
								preSensePair = curSensePair
								agSCAFDebug.debugger.debug("\tprevious sense pair: %s" %(str(preSensePair)))
								agSCAFDebug.debugger.debug("\tcurrent sense pair: %s" %(str(curSensePair)))
							else:
								agSCAFDebug.debugger.debug("\tprevious sense pair: %s" %(str(preSensePair)))
								agSCAFDebug.debugger.debug("\tcurrent sense pair: %s" %(str(curSensePair)))
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
						agSCAFDebug.debugger.debug("\tpossiblePath: %s" %(",".join([vertex2Name[k] for k in possiblePath])))

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
						agSCAFDebug.debugger.debug("\ttmp bestPath: %s" %(",".join([vertex2Name[k] for k in bestPath])))
					agSCAFDebug.debugger.debug("\tBest path: %s" %(",".join([vertex2Name[k] for k in bestPath])))
					scafPaths.append(bestPath)
					agSCAFDebug.debugger.debug("\tbefore %d" %(len(dVisited)))
					for vertex in bestPath:
						visitedLeftOvers.append(k)
					agSCAFDebug.debugger.debug("\tbefore %d" %(len(dVisited)))
					agSCAFDebug.debugger.debug("\tafter %d" %(len(dVisited)))

		agSCAFDebug.debugger.debug("%d visited leftovers %s" %(len(visitedLeftOvers), str([vertex2Name[k] for k in visitedLeftOvers])))

		agSCAFDebug.debugger.debug("number of visited nodes: %d" %(len(dVisited)))
		agSCAFProgress.logger.info("number of visited nodes: %d" %(len(dVisited)))
		agSCAFDebug.debugger.debug("number of paths scaffolded: %d" %(len(scafPaths)))
		return scafPaths

class AGOUTI_Graph(Graph):
	def start(self, joinPairsFile, vertex2Name, dCtgPair2GenePair, algorithm, minSupport):
		self.agSCAFProgress.logger.info("[BEGIN] Scaffolding - Algorithm - %s priority" %(algorithm))
		self.build_graph(joinPairsFile, vertex2Name)
		vertices = self.get_vertices()
		self.agSCAFProgress.logger.info("%d vertices in the graph" %(len(vertices)))

		self.simplify(vertex2Name, minSupport)
		startVertices = self.get_start_vertices()
		self.agSCAFProgress.logger.info("%d vertices to start" %(len(self.get_start_vertices())))

		self.agSCAFProgress.logger.info("Start traversing the GRAPH")
		scafPaths = self.scaffolding(vertex2Name, dCtgPair2GenePair, minSupport)

		self.agSCAFProgress.logger.info("Generating graph in DOT format")
		if self.outGraphFile != "" and self.outGraphFile is not None:
			self.graph2dot(scafPaths, vertex2Name, minSupport)
		self.agSCAFProgress.logger.info("Succeeded")
		return scafPaths

	def get_start_vertices(self):
		self.startVertices = [k for k in self.leaves] + \
							 [k for k in self.nonLeaves if k not in self.leaves]
		return self.startVertices

	def scaffolding(self, vertex2Name, dCtgPair2GenePair, minSupport):
		if self.debug:
			scaffoldingDebug = agLOG.DEBUG("SCAFFOLDING", self.debugLogFile, 'a')

		scafPaths = []
		dVisited = {}

		dSense = {}
		for (vertexA, vertexB), senses in self.senses.iteritems():
			counter = collections.Counter(senses)
			sense = sorted(counter.items(), key=operator.itemgetter(1), reverse=True)[0][0]
#			scaffoldingDebug.debugger.debug("vertexA %s | vertexB %s | %s" %(vertex2Name[vertexA], vertex2Name[vertexB], sense))
			dSense[vertexA, vertexB] = sense

		dBestPaths = {}
		for vertex in self.startVertices:
			if vertex not in dVisited:
				if self.debug:
					scaffoldingDebug.debugger.debug(">%d-%s" %(vertex, vertex2Name[vertex]))
				path = self.dfs(vertex, vertex2Name, minSupport)
				if self.debug:
					scaffoldingDebug.debugger.debug("\t%s" %(",".join([vertex2Name[k] for k in path])))
				vertices = ""
				if len(path) == 1:
					dVisited[vertex] = ""
					if self.debug:
						scaffoldingDebug.debugger.debug("\tpush to as VISITED")
				elif len(path) > 1:
					if self.debug:
						scaffoldingDebug.debugger.debug("\tpush to as VISITED")
					for k in path:
						dVisited[k] = ""
					maxPathLen = 0
					bestPath = []
					for i in range(len(path)):
						if i > 0:
							tmpPath = self.dfs(path[i], vertex2Name, minSupport)
						else:
							tmpPath = path
						scaffoldingDebug.debugger.debug("\ttraverse path from %s: %s"
															 %(vertex2Name[path[i]],
															   ",".join([vertex2Name[k] for k in tmpPath])))

						curVertex = tmpPath[0]
						curCtg = vertex2Name[curVertex]
						curSensePair = ""
						possiblePath = []
						preSensePair = ""
						for i in range(1, len(tmpPath)):
							nextVertex = tmpPath[i]
							nextCtg = vertex2Name[nextVertex]
							if self.debug:
								scaffoldingDebug.debugger.debug("\tcur contig %s | next contig %s"
																%(curCtg, nextCtg))
							if (curVertex, nextVertex) in dSense:
								ctgPair = (curVertex, nextVertex)
								genePair = dCtgPair2GenePair[ctgPair]
								curSensePair = "".join(dSense[ctgPair])
							elif (nextVertex, curVertex) in dSense:
								ctgPair = (nextVertex, curVertex)
								genePair = dCtgPair2GenePair[ctgPair][::-1]
								curSensePair = "".join(dSense[ctgPair][::-1])
								if self.debug:
									scaffoldingDebug.debugger.debug("\treverse contig pair")
							else:
								if self.debug:
									scaffoldingDebug.debugger.debug("\tno connection between this contig pair, stop extension")
								break
#							scaffoldingDebug.debugger.debug("\tgene pair: %s %s"
#															%(genePair[0].geneID, genePair[1].geneID))
							if preSensePair == "":
								preSensePair = curSensePair
								if self.debug:
									scaffoldingDebug.debugger.debug("\tprevious sense pair: %s" %(str(preSensePair)))
									scaffoldingDebug.debugger.debug("\tcurrent sense pair: %s" %(str(curSensePair)))
							else:
								if self.debug:
									scaffoldingDebug.debugger.debug("\tprevious sense pair: %s" %(str(preSensePair)))
									scaffoldingDebug.debugger.debug("\tcurrent sense pair: %s" %(str(curSensePair)))
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
						if self.debug:
							scaffoldingDebug.debugger.debug("\tpossiblePath: %s" %(",".join([vertex2Name[k] for k in possiblePath])))

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
						if self.debug:
							scaffoldingDebug.debugger.debug("\ttmp bestPath: %s"
															%(",".join([vertex2Name[k] for k in bestPath])))
					if self.debug:
						scaffoldingDebug.debugger.debug("\tBest path: %s"
														%(",".join([vertex2Name[k] for k in bestPath])))
					scafPaths.append(bestPath)

		self.agSCAFProgress.logger.info("Number of vertices traversed: %d" %(len(dVisited)))
		self.agSCAFProgress.logger.info("%d paths scaffolded" %(len(scafPaths)))
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

def run_scaffolding(algorithm, vertex2Name, joinPairsFile,
				    dCtgPair2GenePair, outDir, prefix,
					minSupport, debug=0):

	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_scaffolding")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)

	outGraphFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.graph.dot" %(prefix))
	scafPaths = []
	if algorithm == "gene":
		graph = AGOUTI_Graph(outGraphFile)
		graph.start_logger(moduleName, moduleOutDir, prefix, debug)
		scafPaths = graph.start(joinPairsFile, vertex2Name,
								dCtgPair2GenePair, algorithm,
								minSupport)
	elif algorithm == "weight":
		graph = RNAPATHSTAR_Graph(outGraphFile)
		graph.start_logger(moduleName, moduleOutDir, prefix, debug)
		scafPaths = graph.start(joinPairsFile, vertex2Name,
								dCtgPair2GenePair, algorithm,
								minSupport)

	return scafPaths, graph.senses
