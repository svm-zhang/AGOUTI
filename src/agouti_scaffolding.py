import os
import sys
import collections
import itertools
import operator
import time

from lib import agouti_log as agLOG
from src import agouti_path as agPATH

class Graph(object):
	def __init__(self, graph={}):
		self.graph = graph
		self.weights = {}
		self.senses = collections.defaultdict(list)

	def start_logger(self, moduleName, moduleOutDir, prefix, debug=0):
		"""
			initiate a logger to keep track stuff
		"""
		self.debug = debug
		progressLogFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.progressMeter" %(prefix))
		self.agSCAFProgress = agLOG.PROGRESS_METER(moduleName)
		self.agSCAFProgress.add_console_handler()
		self.agSCAFProgress.add_file_handler(progressLogFile)

		self.debugLogFile = None
		if self.debug:
			self.debugLogFile = os.path.join(moduleOutDir, "%s.agouti_scaffolding.debug" %(prefix))

	def build_graph(self, dCtgPairDenoise, vertex2Name):
		"""
			build graph from given joining-pairs
		"""
		self.agSCAFProgress.logger.info("Building graph from joining reads pairs")
		if self.debug:
			buildDebug = agLOG.DEBUG("BUILD", self.debugLogFile)
		for vertexPair, info in dCtgPairDenoise.iteritems():
			vertexA, vertexB = vertexPair
			weight, sense = info
			self.add_vertices(vertexA, vertexB)
			self.add_edge(vertexA, vertexB, weight, sense)

	def add_vertices(self, *vertices):
		"""
			add the give number of contigs
			to the graph
		"""
		for vertex in vertices:
			if vertex not in self.graph:
				self.graph[vertex] = []

	def add_edge(self, vertexA, vertexB, weight, sense):
		"""
			add an edge between the give pair of contig
		"""
		self.graph[vertexA].append(vertexB)
		self.graph[vertexB].append(vertexA)
		self.weights[vertexA, vertexB] = weight
		self.weights[vertexB, vertexA] = weight
		self.senses[vertexA, vertexB] = [sense]*weight

	#def add_edge(self, vertexA, vertexB, senseA, senseB):
	#	"""
	#		add an edge between the give pair of contig
	#	"""
	#	if (vertexA, vertexB) not in self.weights:
	#		self.graph[vertexA].append(vertexB)
	#		#self.weights[vertexA, vertexB] = 0
	#		self.weights[vertexA, vertexB] = 1
	#		#self.senses[vertexA, vertexB].append((senseA, senseB))
	#	else:
	#		self.weights[vertexA, vertexB] += 1
	#	if (vertexB, vertexA) not in self.weights:
	#		self.graph[vertexB].append(vertexA)
	#		#self.weights[vertexB, vertexA] = 0
	#		self.weights[vertexB, vertexA] = 1
	#	else:
	#		self.weights[vertexB, vertexA] += 1
	#	if (vertexA, vertexB) in self.senses:
	#		self.senses[vertexA, vertexB].append((senseA, senseB))
	#	elif (vertexB, vertexA) in self.senses:
	#		self.senses[vertexB, vertexA].append((senseB, senseA))
	#	else:
	#		self.senses[vertexA, vertexB] = [(senseA, senseB)]

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
		# take out each from each other's neighbors
		self.graph[vertexA].remove(vertexB)
		self.graph[vertexB].remove(vertexA)
		# remove the weight between them
		del self.weights[vertexA, vertexB]
		del self.weights[vertexB, vertexA]
		# remove the sense info between them
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
				simplifyDebug.debugger.debug(">vertexA - %d - %s"
											 %(vertexA, vertex2Name[vertexA]))
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
				for i in range(len(edges)):
					weight, vertexB = edges[i]
					if i < 2:
						if weight < minSupport:
							nEdgeRemoved += 1
					else:
						self.remove_edge(vertexA, vertexB)
#				for (weight, vertexB) in edges[3:]:
#					self.remove_edge(vertexA, vertexB)
				if self.debug:
					simplifyDebug.debugger.debug("\tNon-leaf: %s - %s" %(vertexA, vertex2Name[vertexA]))
				self.nonLeaves.append(vertexA)
			elif len(edges) == 1:
				if self.weights[vertexA, edges[0][1]] >= minSupport:
					self.leaves.append(vertexA)
					startVertices.append(vertexA)
					if self.debug:
						simplifyDebug.debugger.debug("\tLeaf: %s - %s" %(vertexA, vertex2Name[vertexA]))
				else:
					nEdgeRemoved += 1
		self.agSCAFProgress.logger.info("%d Edges removed due to insufficient supports"
										%(nEdgeRemoved))

	def graph2dot(self, scafPaths, vertex2Name, minSupport,
				  outDir, prefix):
		"""
			output Graph in DOT format
		"""
		self.agSCAFProgress.logger.info("Visualize graph in DOT")
		outGraphFile = os.path.join(outDir, "%s.agouti_scaffolding.graph.dot" %(prefix))
		with open(outGraphFile, 'w') as fGRAPH:
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

	def report_scaffold_path(self, scafPaths, vertex2Name, outDir, prefix):
		self.agSCAFProgress.logger.info("Report scaffolding paths")
		agPATH.report_scaffold_path(scafPaths, vertex2Name, outDir, prefix)

class AGOUTI_GRAPH_Graph(Graph):
	def start(self, dCtgPairDenoise, vertex2Name, dCtgPair2GenePair, minSupport):
		startTime = time.clock()
		self.build_graph(dCtgPairDenoise, vertex2Name)
		self.agSCAFProgress.logger.info("Build graph took %.4f min CPU time" %((time.clock()-startTime)/60))
		vertices = self.get_vertices()
		self.agSCAFProgress.logger.info("%d vertices in the graph" %(len(vertices)))

		self.simplify(vertex2Name, minSupport)

		self.agSCAFProgress.logger.info("Start graph walk")
		scafPaths = self.scaffolding_v2(vertex2Name, dCtgPair2GenePair, minSupport)
		return scafPaths

	def walk_graph(self, fromVertex, vertex2Name, minSupport, visitedDict={}):
		returnPath = [fromVertex]
		visitedDict[fromVertex] = 1
		toVertex = []
		for toindex in self.graph[fromVertex]:
			weight = self.weights[fromVertex, toindex]
			if weight >= minSupport and toindex not in visitedDict:
				toVertex.append((weight, toindex))

		toVertex.sort(reverse=True)
		for (weight, vertex) in toVertex:
			totalWeight = 0
			for k in self.graph[vertex]:
				totalWeight += self.weights[vertex, k]
			if totalWeight == self.weights[fromVertex,vertex]:
#				self.weights[fromVertex, vertex] = 0
#				self.weights[vertex, fromVertex] = 0
				returnPath += [vertex]
				visitedDict[vertex] = 1
				return returnPath, visitedDict
			else:
#				self.weights[fromVertex, vertex] = 0
#				self.weights[vertex, fromVertex] = 0
				try:
					path, visitedDict = self.walk_graph(vertex, vertex2Name, minSupport, dict(visitedDict, **{vertex:""}))
					returnPath += path
					visitedDict[vertex] = 1
					return returnPath, visitedDict
				except IOError:
					returnPath += [vertex]
					return returnPath, visitedDict
		return returnPath, visitedDict

	def scaffolding_v2(self, vertex2Name, dCtgPair2GenePair, minSupport):
		startTime = time.clock()
		if self.debug:
			scaffoldingDebug = agLOG.DEBUG("SCAFFOLDING", self.debugLogFile, 'a')

		subgraphs = []
		dVisited = {}
		startVertices = self.leaves
		for vertex in startVertices:
			if dVisited.has_key(vertex):
				pass
			else:
				if self.debug:
					scaffoldingDebug.debugger.debug(">graph walk start from node: %s"
													%(vertex2Name[vertex]))
				subgraph, dVisited = self.walk_graph(vertex, vertex2Name, minSupport, dVisited) #added
				if len(subgraph) > 1:
					if self.debug:
						scaffoldingDebug.debugger.debug("return subgraph: %s"
														%(",".join(map(str, [vertex2Name[k] for k in subgraph]))))
						scaffoldingDebug.debugger.debug("number of nodes in this subgraph: %d"
														%(len(subgraph)))
					subgraphs.append(subgraph)
		if self.debug:
			scaffoldingDebug.debugger.debug("number of subgraphs: %d"
											%(len(subgraphs)))

		orphans = [k for k in self.get_vertices() if k not in dVisited]
		orphans += [k for k in self.nonLeaves if k not in dVisited and k not in orphans]
		if self.debug:
			scaffoldingDebug.debugger.debug("Number of orphan vertices left: %d " %(len(orphans)))
			scaffoldingDebug.debugger.debug("Finding paths from non-leaves")

		for vertex in orphans:
			if vertex not in dVisited:
				subgraph, dVisited = self.walk_graph(vertex, vertex2Name, minSupport, dVisited)
				if self.debug:
					scaffoldingDebug.debugger.debug(">graph walk start from node: %s"
													%(vertex2Name[vertex]))
					scaffoldingDebug.debugger.debug("missing subgraph: %s"
													%(",".join([vertex2Name[k] for k in subgraph])))
#				vertices = ""
#				for k in subgraph:
#					if k in dVisited:
#						scaffoldingDebug.debugger.debug("seen %s before" %(vertex2Name[k]))
#						subgraph.remove(k)
#					else:
#						dVisited[k] = 1
				if len(subgraph) > 1:
					if self.debug:
						scaffoldingDebug.debugger.debug("add to subgraphs: %s"
														%(str([vertex2Name[k] for k in subgraph])))
					subgraphs.append(subgraph)
		if self.debug:
			scaffoldingDebug.debugger.debug("number of visited nodes: %d" %(len(dVisited)))
		self.agSCAFProgress.logger.info("number of visited nodes: %d" %(len(dVisited)))
		self.agSCAFProgress.logger.info("Scaffolding took %.4f min CPU time" %((time.clock()-startTime)/60))

		dSense = {}
		for (fVertex, tVertex), senses in self.senses.iteritems():
			counter = collections.Counter(senses)
			sense = sorted(counter.items(), key=operator.itemgetter(1), reverse=True)[0][0]
			dSense[fVertex, tVertex] = sense

		startTime = time.clock()
		if self.debug:
			scaffoldingDebug.debugger.debug("Graph Reconciliation")
		self.agSCAFProgress.logger.info("Graph Reconciliation")
		scafPaths = []
		for path in subgraphs:
			if self.debug:
				scaffoldingDebug.debugger.debug(">subGraphs: %s" %(str([vertex2Name[k] for k in path])))
				bestOrder = self.reconcile(path, dVisited, dCtgPair2GenePair,
										   dSense, vertex2Name, minSupport,
										   scaffoldingDebug)
			else:
				bestOrder = self.reconcile(path, dVisited, dCtgPair2GenePair,
										   dSense, vertex2Name, minSupport)
			if len(bestOrder) > 1:
				scafPaths.append(bestOrder)
		self.agSCAFProgress.logger.info("Reconciliation took %.4f min CPU time" %((time.clock()-startTime)/60))

		if self.debug:
			scaffoldingDebug.debugger.debug("number of paths scaffolded: %d" %(len(scafPaths)))
		return scafPaths

	def dfs(self, vertex, vertex2Name, minSupport, seen=None, path=None):
		"""
			Depth-first traversal from a given node
		"""
		if seen is None: seen = [vertex]
		if path is None: path = [vertex]

		neighbors = []
		for vertexB in self.graph[vertex]:
			weight = self.weights[vertex, vertexB]
			neighbors.append((weight, vertexB))

		neighbors.sort(reverse=True)
		for weight, vertexB in neighbors:
			if weight >= minSupport and vertexB not in seen:
				seen.append(vertexB)
				path.append(vertexB)
				path = self.dfs(vertexB, vertex2Name, minSupport, seen, path)
		return path

	def reconcile(self, path, dVisited, dCtgPair2GenePair,
				  dSense, vertex2Name, minSupport,
				  scaffoldingDebug=None):
		maxPathLen = 0
		dWeight = {}
		bestPath = []
		for i in range(len(path)):
			if i > 0:
				tmpPath, _ = self.walk_graph(path[i], vertex2Name, minSupport)
			else:
				tmpPath = path
			if self.debug:
				scaffoldingDebug.debugger.debug("\twalking path: %s" %(",".join([vertex2Name[k] for k in tmpPath])))
			curVertex = tmpPath[0]
			curCtg = vertex2Name[curVertex]
			totalWeight = 0
			curSensePair = ""
			possiblePath = []
			preSensePair = ""
			for i in range(1, len(tmpPath)):
				nextVertex = tmpPath[i]
				nextCtg = vertex2Name[nextVertex]
				if self.debug:
					scaffoldingDebug.debugger.debug("\t\tcur %s | next %s" %(curCtg, nextCtg))
				if (curVertex, nextVertex) in self.weights:
					if (curVertex, nextVertex) in dSense:
						ctgPair = (curVertex, nextVertex)
						genePair = dCtgPair2GenePair[ctgPair]
						curSensePair = "".join(dSense[ctgPair])
						weight = self.weights[ctgPair]
					elif (nextVertex, curVertex) in dSense:
						ctgPair = (nextVertex, curVertex)
						genePair = dCtgPair2GenePair[ctgPair][::-1]
						curSensePair = "".join(dSense[ctgPair][::-1])
						if self.debug:
							scaffoldingDebug.debugger.debug("\t\treverse contig pair")
						weight = self.weights[ctgPair]
					else:
						if self.debug:
							scaffoldingDebug.debugger.debug("\t\tthe pair has been removed, stop extension")
						break
				else:
					if self.debug:
						scaffoldingDebug.debugger.debug("\t\tthe edge between these two contigs do not exist")
					break
				if self.debug:
					scaffoldingDebug.debugger.debug("\t\tgene pair: %s %s" %(genePair[0].geneID, genePair[1].geneID))
				if preSensePair == "":
					if self.debug:
						scaffoldingDebug.debugger.debug("\t\tprevious sense pair: %s" %(str(preSensePair)))
						scaffoldingDebug.debugger.debug("\t\tcurrent sense pair: %s" %(str(curSensePair)))
				else:
					if self.debug:
						scaffoldingDebug.debugger.debug("\t\tprevious sense pair: %s" %(str(preSensePair)))
						scaffoldingDebug.debugger.debug("\t\tcurrent sense pair: %s" %(str(curSensePair)))
					if preSensePair == "+-" and (curSensePair == "+-" or curSensePair == "++"):
						pass
					elif preSensePair == "++" and (curSensePair == "--" or curSensePair == "-+"):
						pass
					elif preSensePair == "--" and (curSensePair == "+-" or curSensePair == "++"):
						pass
					elif preSensePair == "-+" and (curSensePair == "-+" or curSensePair == "--"):
						pass
					else:
						if self.debug:
							scaffoldingDebug.debugger.debug("\t\tExtension terminated")
						break
				totalWeight += weight
				possiblePath.append(curVertex)
				preSensePair = curSensePair
				curVertex = nextVertex
				curCtg = vertex2Name[curVertex]
			possiblePath.append(curVertex)
			if self.debug:
				scaffoldingDebug.debugger.debug("\tpossiblePath: %s"
												%(",".join([vertex2Name[k] for k in possiblePath])))
			if len(possiblePath) == len(path):
				bestPath = possiblePath
				break
			else:
				if maxPathLen == 0:
					maxPathLen = len(possiblePath)
					bestPath = possiblePath
					maxWeight = totalWeight
				else:
					if maxPathLen < len(possiblePath):
						maxPathLen = len(possiblePath)
						bestPath = possiblePath
					elif maxPathLen == len(possiblePath):
						if maxWeight < totalWeight:
							bestPath = possiblePath
			if self.debug:
				scaffoldingDebug.debugger.debug("\ttmp bestPath: %s" %(",".join([vertex2Name[k] for k in bestPath])))
		if self.debug:
			scaffoldingDebug.debugger.debug("\tBest path: %s" %(",".join([vertex2Name[k] for k in bestPath])))
		return bestPath

def run_scaffolding(vertex2Name, dCtgPairDenoise,
				    dCtgPair2GenePair, outDir, prefix,
					minSupport, debug=0):

	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_scaffolding")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)

#	if algorithm == "gene":
#		graph = AGOUTI_Graph(outGraphFile)
#		graph.start_logger(moduleName, moduleOutDir, prefix, debug)
#		scafPaths = graph.start(joinPairsFile, vertex2Name,
#								dCtgPair2GenePair, algorithm,
#								minSupport)
	graph = AGOUTI_GRAPH_Graph()
	graph.start_logger(moduleName, moduleOutDir, prefix, debug)
	scafPaths = graph.start(dCtgPairDenoise, vertex2Name,
							dCtgPair2GenePair,
							minSupport)
	graph.report_scaffold_path(scafPaths, vertex2Name, outDir, prefix)
	graph.graph2dot(scafPaths, vertex2Name, minSupport,
					outDir, prefix)

	return scafPaths, graph.senses
