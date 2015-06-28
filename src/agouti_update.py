import os
import sys
import collections

from lib import agouti_gff as agGFF

def agouti_update(pathList, contigDict, nameList, origSize,
				  edgeSenseDict, visitedDict, dGFFs,
				  dCtgPair2GenePair, outDir, prefix,
				  lenGapFill=1000):
#				  outFasta, outScaffPath, lenGapFill=1000):

	sys.stderr.write("Updating gene models ... \n")
	scafID = 0
	newSizeList = []
	scafName = ""
	outFasta = os.path.join(outDir, "%s.agouti.fasta" %(prefix))
	outGFF = os.path.join(outDir, "%s.agouti.gff" %(prefix))
	outScaffPath = os.path.join(outDir, "%s.scaff.paths" %(prefix))
	fOUTFASTA = open(outFasta, "w")
	fSCAFFPATH = open(outScaffPath, "w")
	dUpdateGFFs = collections.defaultdict(list)
	numMergedGene = 0
	numScaffolded = 0
	numSkipped = 0		# pairs of contigs in the path do not exist in dCtgPair2GenePair

#	flag = 0
#	dProblem = []
#	print "updating ...\n"
#	for i in range(len(pathList)):
#		path = pathList[i]
#		curNode = path[0]
#		for node in path[1:]:
#			if node - curNode > 2:
#				flag = 1
#				break
#			curNode = node
#
#		if flag == 1:
#			dProblem.append(i)
#			flag = 0

#	numVisitedNodes = sum([len(path) for path in pathList])
#	print "numVisitedNodes", numVisitedNodes
#	print "visitedDict", len(visitedDict)
#	twice = {}
#	for path in pathList:
#		for node in path:
#			if node in visitedDict:
#				if node not in twice:
#					twice[node] = 1
#				else:
#					twice[node] += 1
#	for k, v in twice.iteritems():
#		if v == 2:
#			print k, nameList[k]

	loop2Remove, visitedDict = remove_cycles(pathList, visitedDict)

	print "updating ... "
	for i in range(len(pathList)):
		path = pathList[i]
		if len(path) >= 3 and i not in loop2Remove:
			tmp = 0
			scafID += 1
			scafName = prefix + "_scaf_%d" %(scafID)
			fSCAFFPATH.write("%s: %s\n" % (scafName, str(path)))
			vertexNameList = []
			for vertex in path:
				vertexNameList.append(nameList[vertex])
			pathDescription = ",".join(vertexNameList)

			print >> fSCAFFPATH, pathDescription
			currentVertex = path[0]
			sequence = contigDict[currentVertex]
			currentSense = "+"
			currentCtg = nameList[currentVertex]
			assemblyList = currentVertex
			genePath = ()
			preGeneID, curGeneID = "", ""
			mergedGene = agGFF.AGOUTI_GFF()
			preMergedGene = agGFF.AGOUTI_GFF()
			gapStart, gapStop = 0, 0
			offset = 0
			orientation = ""
			updatedGeneIDs = []
			excludeGeneIDs = []
			for nextVertex in path[1:]:
				nextCtg = nameList[nextVertex]

				if tmp == 0:
					print "scaf_%d" %(scafID), path
					tmp = 1
				print ">>", currentVertex, nextVertex, currentCtg, nextCtg

				currentGene, nextGene = ctgpair2genepair(dCtgPair2GenePair, currentCtg, nextCtg)
				if currentGene is None and nextGene is None:
					print "%s %s not found in dCtgPair2GenePair" %(currentCtg, nextCtg)
					#!!! I should not break here, should continue#
					numSkipped += 1
					break
				curGeneID = currentGene.geneID
				excludeGeneIDs = [preGeneID] + [curGeneID]
				print ">>>>", scafName, "pre", preGeneID, "current", currentGene.geneID, "next", nextGene.geneID

				toLeft, toRight = 0, 0
				if preGeneID != "":
					if curGeneID != preGeneID:
						toLeft, toRight = getGeneIndex(dGFFs[currentCtg]+dUpdateGFFs[scafName], curGeneID)
					else:
						# problematic
#						if len(dGFFs[currentCtg]) == 1:
						toLeft, toRight = getGeneIndex(dUpdateGFFs[scafName], preMergedGene.geneID)
#						else:
#							toLeft, toRight = getGeneIndex(dGFFs[currentCtg]+dUpdateGFFs[scafName], currentGene.geneID)
				print ">>>>", "toleft", toLeft, "toRight", toRight

				FR, FF, RR, RF = get_orientation_counts(currentVertex, nextVertex, edgeSenseDict)
				print ">>>>orientation", currentVertex, nextVertex, currentCtg, nextCtg, currentSense, FF, RR, RF, FR
				if currentSense == "-":
					# we had flipped the upstream piece! Must flip again
					temp1 = FR
					temp2 = FF
					FR = RR
					FF = RF
					RR = temp1
					RF = temp2
				orientation = get_scaffolding_orientation(FR, FF, RR, RF)

				print ">>>>orientation", currentVertex, nextVertex, currentCtg, nextCtg, currentSense, FF, RR, RF, FR
				gapStart = gapStop + len(contigDict[currentVertex]) + 1
				gapStop = gapStart + lenGapFill - 1
				print ">>>>", "offset", offset, len(contigDict[currentVertex]), "gapstart", gapStart, "gapstop", gapStop+1
				if orientation == "FR":
					if curGeneID !=  preGeneID:
						numMergedGene += 1
						mergedGene = merge_gene_model(currentGene, nextGene, scafName,
													  numMergedGene, offset, gapStop+1)
						dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
					else:
						mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
													  numMergedGene, 0, gapStop+1)
						indexMerged = updatedGeneIDs.index(mergedGene.geneID)
						dUpdateGFFs[scafName][indexMerged] = mergedGene
					preMergedGene = mergedGene
					assemblyList = ((assemblyList, "+"), (nextVertex, "+"))
					sequence += 'N'*lenGapFill + contigDict[nextVertex]
					currentSense = "+"
				elif orientation == "FF":
					nextGene = reverse_gene_model(nextGene, len(contigDict[nextVertex]))
					if curGeneID !=  preGeneID:
						numMergedGene += 1
						mergedGene = merge_gene_model(currentGene, nextGene, scafName,
													  numMergedGene, offset, gapStop+1)
						dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
					else:
						mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
													  numMergedGene, 0, gapStop+1)
						indexMerged = updatedGeneIDs.index(mergedGene.geneID)
						dUpdateGFFs[scafName][indexMerged] = mergedGene
					preMergedGene = mergedGene
					assemblyList = ((assemblyList, "+"), (nextVertex, "-"))
					sequence += 'N'*lenGapFill + complement(contigDict[nextVertex])
					currentSense = "-"
				elif orientation == "RR":
					if toLeft <= 1:
						dGFFs[currentCtg] = reverse_gene_models(dGFFs[currentCtg], len(contigDict[currentVertex]))
						if currentGene.geneID != preGeneID:
							numMergedGene += 1
							mergedGene = merge_gene_model(currentGene, nextGene, scafName,
														  numMergedGene, offset, gapStop+1)
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
						else:
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs)
							dUpdateGFFs[scafName] = reverse_gene_models(dUpdateGFFs[scafName], gapStart)
							mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
														  numMergedGene, 0, gapStop+1)
							indexMerged = updatedGeneIDs.index(mergedGene.geneID)
							dUpdateGFFs[scafName][indexMerged] = mergedGene
						assemblyList = ((assemblyList, "-"), (nextVertex, "+"))
						sequence = complement(contigDict[currentVertex]) + \
								   'N'*lenGapFill + contigDict[nextVertex]
						preMergedGene = mergedGene
						print ">>>> RR zone preMergedGene", preMergedGene.geneID, preMergedGene.geneStart, preMergedGene.geneStop
					else:
						currentVertex = nextVertex
						if path.index(currentVertex) < len(path)-1:
							sequence = contigDict[currentVertex]
							currentCtg = nextCtg
							currentSense = "+"
							assemblyList = currentVertex
							scafID += 1
							scafName = prefix + "_scaf_%d" %(scafID)
							gapStart, gapStop = 0, 0
							offset = 0
							orientation = ""
							updatedGeneIDs = []
							excludeGeneIDs = []
							preGeneID = ""
						else:
							break
					currentSense = "+"
				elif orientation == "RF":
					if toLeft <= 1:
						if currentGene.geneID != preGeneID:
							dGFFs[currentCtg] = reverse_gene_models(dGFFs[currentCtg], len(contigDict[currentVertex]))
							nextGene = reverse_gene_model(nextGene, len(contigDict[nextVertex]))
							numMergedGene += 1
							mergedGene = merge_gene_model(currentGene, nextGene, scafName,
														  numMergedGene, offset, gapStop+1)
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
						else:
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs)
							dUpdateGFFs[scafName] = reverse_gene_models(dUpdateGFFs[scafName], gapStop+len(contigDict[currentVertex]))
							mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
														  numMergedGene, 0, gapStop+1)
							indexMerged = updatedGeneIDs.index(mergedGene.geneID)
							dUpdateGFFs[scafName][indexMerged] = mergedGene
						assemblyList = ((assemblyList, "-"), (nextVertex, "-"))
						sequence = complement(contigDict[currentVertex]) + \
								   'N'*lenGapFill + \
								   complement(contigDict[nextVertex])
						preMergedGene = mergedGene
						print ">>>> RF zone preMergedGene", preMergedGene.geneID, preMergedGene.geneStart, preMergedGene.geneStop
					else:
						currentVertex = nextVertex
						if path.index(currentVertex) < len(path)-1:
							sequence = contigDict[currentVertex]
							currentCtg = nextCtg
							currentSense = "+"
							assemblyList = currentVertex
							scafID += 1
							scafName = prefix + "_scaf_%d" %(scafID)
							gapStart, gapStop = 0, 0
							offset = 0
							orientation = ""
							updatedGeneIDs = []
							excludeGeneIDs = []
							preGeneID = ""
						else:
							break
					currentSense = "-"
				offset = gapStop
				preGeneID = nextGene.geneID
				currentVertex = nextVertex
				currentCtg = nameList[currentVertex]

				outstring = ">>>>(%d, %d): FF %d RR %d RF %d FR %d : %s \t%s" % (currentVertex, nextVertex, FF, RR, RF, FR, orientation, str(assemblyList))
				print outstring
				print >> fSCAFFPATH, outstring

			excludeGeneIDs = [preGeneID]
			dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs)
			fOUTFASTA.write(">%s|%dbp|%s|%s\n%s\n" % (scafName, len(sequence), pathDescription, str(assemblyList), sequence))
			print
	fSCAFFPATH.close()

	# other contigs need to be output
	sys.stderr.write("outputting contigs spared from scaffolding ...\n")
	numLeft = 0
	for vertex in contigDict:
		if vertex in visitedDict:
				if nameList[vertex] in dGFFs:
					del dGFFs[nameList[vertex]]
				continue
		numLeft += 1
		fOUTFASTA.write(">%s\n%s\n" % (nameList[vertex], contigDict[vertex]))
	fOUTFASTA.close()

	# output in GFF format
	sys.stderr.write("outputting Gene Models in GFF3 ...\n")
	dFinalGFFs = dict(dGFFs, **dUpdateGFFs)
	numGenes = gff_outputter(dFinalGFFs, outGFF)

	sys.stderr.write("\n####Summary####\n")

	sys.stderr.write("number of contigs scaffoled: %d\n" %(len(visitedDict)))
	sys.stderr.write("number of scaffolds: %d\n" %(scafID))
	sys.stderr.write("number of contigs found no links: %d\n" %(numLeft))
	sys.stderr.write("number of contigs in the final assembly: %d\n" %(numLeft+scafID))
	sys.stderr.write("Final number of genes: %d\n" %(numGenes))

def gff_outputter(dGeneModels, outGFF):
	fOUTGFF = open(outGFF, 'w')
	fOUTGFF.write("##gff-version3\n")
	fOUTGFF.write("# This output was generated with AGOUTI (version 0.1)\n")
	numGenes = 0
	for k, v in dGeneModels.iteritems():
		numGenes += len(v)
		for i in range(len(v)):
			geneModel = v[i]
			fOUTGFF.write("# start gene %s\n" %(geneModel.geneID))
			fOUTGFF.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s\n" %(k, geneModel.program, geneModel.geneStart, geneModel.geneStop, geneModel.strand, geneModel.geneID))
			for j in range(0, len(geneModel.lcds), 2):
				cdsStart = geneModel.lcds[j]
				cdsStop = geneModel.lcds[j+1]
				fOUTGFF.write("%s\t%s\tCDS\t%d\t%d\t.\t%s\t.\tID=%s.cds\n" %(k, geneModel.program, cdsStart, cdsStop, geneModel.strand, geneModel.geneID))
			fOUTGFF.write("# end gene %s\n" %(geneModel.geneID))
			fOUTGFF.write("###\n")
	fOUTGFF.close()
	return numGenes

def remove_cycles(pathList, visitedDict):
	'''
		removing cyclic paths
		Still need to figure out why this can happen
	'''
	loop2Remove = []
	loop = 0
	for i in range(len(pathList)):
		path = pathList[i]
		nodes = []
		for node in path:
			if node not in nodes:
				nodes.append(node)
			else:
				loop = 1
				break
		if loop == 1:
			loop2Remove.append(i)
			loop = 0
			for node in path:
				if node in visitedDict:
					del visitedDict[node]
	return loop2Remove, visitedDict

def getGeneIndex(geneModels, curGeneID):
	currentGeneIndex = [x.geneID for x in geneModels].index(curGeneID)
	toLeft = currentGeneIndex
	toRight = len(geneModels) - currentGeneIndex
	return toLeft, toRight

def merge_gene_model(currentGene, nextGene, scafName,
					 numMergedGene, startOffset, stopOffset):
	mergedGene = agGFF.AGOUTI_GFF()
	mergedGene.geneID = "AGOUTI_Merged_gene_%d" %(numMergedGene)
	mergedGene.geneStart = currentGene.geneStart + startOffset
	mergedGene.geneStop = nextGene.geneStop + stopOffset
	mergedGene.program = "AGOUTI"
	mergedGene.ctgID = scafName
	if currentGene.strand == nextGene.strand:
		mergedGene.strand = currentGene.strand
	else:
		mergedGene.strand = '.'
	print ">>>> merge", "current LCDS", currentGene.lcds
	for i in range(len(currentGene.lcds)):
		currentGene.lcds[i] += startOffset
	print ">>>> merge", "next LCDS", nextGene.lcds
	for i in range(len(nextGene.lcds)):
		nextGene.lcds[i] += stopOffset
	mergedGene.lcds = currentGene.lcds + nextGene.lcds
	return mergedGene

def update_gene_model(geneModels, updatedGeneModels,
					  newScafID, offset, excludes, 
					  mergedGene=None):
	print ">>>>", excludes
	for i in range(len(geneModels)):
		geneModels[i].ctgID = newScafID
		geneModels[i].geneStart += offset
		geneModels[i].geneStop += offset
		if geneModels[i].geneID not in excludes:
			updatedGeneModels.append(geneModels[i])
			for j in range(len(geneModels[i].lcds)):
				geneModels[i].lcds[j] += offset
	if mergedGene is not None:
		updatedGeneModels.append(mergedGene)
	updatedGeneIDs = []
	for i in range(len(updatedGeneModels)):
		updatedGeneIDs.append(updatedGeneModels[i].geneID)
	return updatedGeneModels, updatedGeneIDs

def reverse_gene_models(geneModels, lenCtg, excludeIDs=[]):
	print ">>>>reverse, %d" %(lenCtg)
	for i in xrange(len(geneModels)):
#		geneModel = geneModels[i]
		if geneModels[i].geneID not in excludeIDs:
			geneModels[i] = reverse_gene_model(geneModels[i], lenCtg)
	return geneModels

def reverse_gene_model(geneModel, lenCtg):
	tmpGeneStop =  lenCtg-geneModel.geneStart+1
	tmpGeneStart = lenCtg-geneModel.geneStop+1
	geneModel.geneStart = tmpGeneStart
	geneModel.geneStop = tmpGeneStop
	lcds = geneModel.lcds
	reverseLCDs = []
	print ">>>> before reverse", "LCDS", lcds
	for j in range(len(lcds)-1, -1, -2):
		reverseLCDs += [lenCtg-lcds[j]+1, lenCtg-lcds[j-1]+1]
	geneModel.lcds = reverseLCDs
	if geneModel.strand == '+':
		geneModel.strand = '-'
	else:
		geneModel.strand = '+'
#	print lcds
#	print reverseLCDs
	return geneModel

def ctgpair2genepair(dCtgPair2GenePair, currentCtg, nextCtg):
	currentGene = None
	nextGene = None
	if (currentCtg, nextCtg) in dCtgPair2GenePair:
		genePair = dCtgPair2GenePair[currentCtg, nextCtg]
		currentGene = genePair[0]
		nextGene = genePair[1]
	elif (nextCtg, currentCtg) in dCtgPair2GenePair:
		genePair = dCtgPair2GenePair[nextCtg, currentCtg]
		currentGene = genePair[1]
		nextGene = genePair[0]
	return currentGene, nextGene

def get_orientation_counts(currentVertex, nextVertex, edgeSenseDict):
	''' get the counts of each orientation '''
	senseList = []
	if (currentVertex, nextVertex) in edgeSenseDict:
		senseList = edgeSenseDict[currentVertex, nextVertex]
		FR = senseList.count(("+", "-"))
		RF = senseList.count(("-", "+"))
	elif (nextVertex, currentVertex) in edgeSenseDict:
		senseList = edgeSenseDict[nextVertex, currentVertex]
		# flip because the from and end vertices fliped
		FR = senseList.count(("-", "+"))
		RF = senseList.count(("+", "-"))
	else:
		print "do not get strand orientation for these two nodes: %d, %d" %(currentVertex, nextVertex)
		sys.exit()
	FF = senseList.count(("+", "+"))
	RR = senseList.count(("-", "-"))
	return FR, FF, RR, RF

def get_scaffolding_orientation(FR, FF, RR, RF):
	''' get scaffolding orientation '''

	# we have FR - leave alone
	if FR >= FF and FR >= RR and FR >= RF:
		return "FR"
	# we have FF - flip righthand side
	elif FF >= RR and FF >= RF:
		return "FF"
	# we have RR - flip lefthand side
	elif RR >= RF:
		return "RR"
	else:
		return "RF"

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
