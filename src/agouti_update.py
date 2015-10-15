import os
import sys
import collections

from lib import agouti_gff as agGFF

def agouti_update(pathList, contigDict, nameList,
				  edgeSenseDict, visitedDict, dGFFs,
				  dCtgPair2GenePair, outDir, prefix, numNs=1000):
	sys.stderr.write("Updating gene models ... \n")
	scafID = 0

	print "updating ... "
	outScaffPath = os.path.join(outDir, "%s.scaff.paths" %(prefix))
	fSCAFFPATH = open(outScaffPath, "w")
	outFasta = os.path.join(outDir, "%s.agouti.fasta" %(prefix))
	fOUTFASTA = open(outFasta, "w")
	dUpdateGFFs = collections.defaultdict(list)
	dMergedGene2Ctgs = collections.defaultdict(list)
	dMergedGene2Genes = collections.defaultdict(list)
	numMergedGene = 0
	nCtgScaffolded = 0
	for i in range(len(pathList)):
		path = pathList[i]
		if len(path) >= 2:
			tmp = 0
			scafID += 1
			scafName = prefix + "_scaf_%d" %(scafID)
			fSCAFFPATH.write(">%s\n" % (scafName))

			currentVertex = path[0]
			sequence = contigDict[currentVertex]
			currentSense = "+"
			currentCtg = nameList[currentVertex]
			scafPath = []
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
			stopFlag = 0
			for nextVertex in path[1:]:
				print "sequence length", len(sequence)
				nextCtg = nameList[nextVertex]

				if tmp == 0:
					print "scaf_%d" %(scafID), path
					tmp = 1
				print ">>", currentVertex, nextVertex, currentCtg, nextCtg

				currentGene, nextGene = ctgpair2genepair(dCtgPair2GenePair, currentCtg, nextCtg)
				if currentGene is None and nextGene is None:
					print "%s %s not found in dCtgPair2GenePair" %(currentCtg, nextCtg)
					#!!! I should not break here, should continue#
					break
				curGeneID = currentGene.geneID
				excludeGeneIDs = [preGeneID] + [curGeneID]
				print ">>>>", scafName, "pre", preGeneID, "current", currentGene.geneID, "next", nextGene.geneID

				toLeft, toRight = 0, 0
				if preGeneID != "":
					if curGeneID != preGeneID:
						if currentSense == "-":
							tmp = exclude_gene_model(dGFFs[currentCtg], preGeneID)[::-1]
							print "tmp", " ".join([x.geneID for x in tmp])
							toLeft, toRight = get_gene_index(dUpdateGFFs[scafName]+tmp, curGeneID)
						else:
							toLeft, toRight = get_gene_index(dUpdateGFFs[scafName]+dGFFs[currentCtg], curGeneID)
					else:
						# problematic
#						if len(dGFFs[currentCtg]) == 1:
						tmp = exclude_gene_model(dGFFs[currentCtg], preGeneID)
						print "tmp", " ".join([x.geneID for x in tmp])
						toLeft, toRight = get_gene_index(dUpdateGFFs[scafName]+tmp, preMergedGene.geneID)
#						else:
#							toLeft, toRight = get_gene_index(dGFFs[currentCtg]+dUpdateGFFs[scafName], currentGene.geneID)
				print ">>>>", "toleft", toLeft, "toRight", toRight

				FR, FF, RR, RF = get_orientation_counts(currentVertex, nextVertex, edgeSenseDict)
				print ">>>>orientation", currentVertex, nextVertex, currentCtg, nextCtg, currentSense, "FF", FF, "RR", RR, "RF", RF, "FR", FR
				if currentSense == "-":
					# we had flipped the upstream piece! Must flip again
					temp1 = FR
					temp2 = FF
					FR = RR
					FF = RF
					RR = temp1
					RF = temp2
				orientation = get_scaffolding_orientation(FR, FF, RR, RF)

				print ">>>>orientation", currentVertex, nextVertex, currentCtg, nextCtg, currentSense, "FF", FF, "RR", RR, "RF", RF, "FR", FR
				gapStart = gapStop + 1 + len(contigDict[currentVertex])
				gapStop = gapStart + numNs - 1
				print ">>>>", "offset", offset, "current Contig length", len(contigDict[currentVertex]), "gapstart", gapStart, "gapstop", gapStop+1
				if orientation == "FR":
					if curGeneID !=  preGeneID:
						numMergedGene += 1
						mergedGene = merge_gene_model(currentGene, nextGene, scafName,
													  numMergedGene, offset, gapStop)
						dMergedGene2Ctgs[mergedGene.geneID] += [currentCtg, nextCtg]
						dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
						dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
					else:
						mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
													  numMergedGene, 0, gapStop)
						dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
						dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
						indexMerged = updatedGeneIDs.index(mergedGene.geneID)
						dUpdateGFFs[scafName][indexMerged] = mergedGene
					print "mergedGene", mergedGene.geneID, "dMergedGene2Ctgs", dMergedGene2Ctgs[mergedGene.geneID]
					preMergedGene = mergedGene
					assemblyList = ((assemblyList, "+"), (nextVertex, "+"))
					sequence += 'N'*numNs + contigDict[nextVertex]
					currentSense = "+"
				elif orientation == "FF":
					nextGene = reverse_gene_model(nextGene, len(contigDict[nextVertex]))
					if curGeneID !=  preGeneID:
						numMergedGene += 1
						mergedGene = merge_gene_model(currentGene, nextGene, scafName,
													  numMergedGene, offset, gapStop)
						dMergedGene2Ctgs[mergedGene.geneID] += [currentCtg, nextCtg]
						dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
						dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
					else:
						mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
													  numMergedGene, 0, gapStop)
						dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
						dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
						indexMerged = updatedGeneIDs.index(mergedGene.geneID)
						dUpdateGFFs[scafName][indexMerged] = mergedGene
					preMergedGene = mergedGene
					assemblyList = ((assemblyList, "+"), (nextVertex, "-"))
					sequence += 'N'*numNs + rc_seq(contigDict[nextVertex])
					currentSense = "-"
				elif orientation == "RR":
					if toLeft < 1:
						if currentGene.geneID != preGeneID:
							dGFFs[currentCtg] = reverse_gene_models(dGFFs[currentCtg], len(contigDict[currentVertex]))
							numMergedGene += 1
							mergedGene = merge_gene_model(currentGene, nextGene, scafName,
														  numMergedGene, offset, gapStop)
							dMergedGene2Ctgs[mergedGene.geneID] += [currentCtg, nextCtg]
							dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
						else:
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs)
							dUpdateGFFs[scafName] = reverse_gene_models(dUpdateGFFs[scafName], gapStart-1)
							mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
														  numMergedGene, 0, gapStop)
							dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
							dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
							indexMerged = updatedGeneIDs.index(mergedGene.geneID)
							dUpdateGFFs[scafName][indexMerged] = mergedGene
						assemblyList = ((assemblyList, "-"), (nextVertex, "+"))
						sequence = rc_seq(sequence) + \
								   'N'*numNs + contigDict[nextVertex]
						preMergedGene = mergedGene
						print ">>>> RR zone preMergedGene", preMergedGene.geneID, preMergedGene.geneStart, preMergedGene.geneStop
					else:
						currentVertex = nextVertex
						if path.index(currentVertex) < len(path)-1:
							sequence = contigDict[currentVertex]
							currentCtg = nextCtg
							currentSense = "+"
							assemblyList = currentVertex
							print "delete before", " ".join([x.geneID for x in dUpdateGFFs[scafName]])
							del dUpdateGFFs[scafName]
							scafPath = []
							scafName = prefix + "_scaf_%d" %(scafID)
							gapStart, gapStop = 0, 0
							offset = 0
							orientation = ""
							updatedGeneIDs = []
							excludeGeneIDs = []
							preGeneID = ""
							continue
						else:
							print "RR", "here is a stop"
							break
					currentSense = "+"
				elif orientation == "RF":
					if toLeft < 1:
						if currentGene.geneID != preGeneID:
							dGFFs[currentCtg] = reverse_gene_models(dGFFs[currentCtg], len(contigDict[currentVertex]))
							nextGene = reverse_gene_model(nextGene, len(contigDict[nextVertex]))
							numMergedGene += 1
							mergedGene = merge_gene_model(currentGene, nextGene, scafName,
														  numMergedGene, offset, gapStop)
							dMergedGene2Ctgs[mergedGene.geneID] += [currentCtg, nextCtg]
							dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs, mergedGene)
						else:
							dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs)
							dUpdateGFFs[scafName] = reverse_gene_models(dUpdateGFFs[scafName], gapStop+len(contigDict[currentVertex]))
							mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
														  numMergedGene, 0, gapStop)
							dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
							dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
							indexMerged = updatedGeneIDs.index(mergedGene.geneID)
							dUpdateGFFs[scafName][indexMerged] = mergedGene
						assemblyList = ((assemblyList, "-"), (nextVertex, "-"))
						sequence = rc_seq(sequence) + \
								   'N'*numNs + \
								   rc_seq(contigDict[nextVertex])
						preMergedGene = mergedGene
						print ">>>> RF zone preMergedGene", preMergedGene.geneID, preMergedGene.geneStart, preMergedGene.geneStop
					else:
						currentVertex = nextVertex
						if path.index(currentVertex) < len(path)-1:
							sequence = contigDict[currentVertex]
							currentCtg = nextCtg
							currentSense = "+"
							assemblyList = currentVertex
							print "delete before", " ".join([x.geneID for x in dUpdateGFFs[scafName]])
							del dUpdateGFFs[scafName]
							scafPath = []
							scafName = prefix + "_scaf_%d" %(scafID)
							gapStart, gapStop = 0, 0
							offset = 0
							orientation = ""
							updatedGeneIDs = []
							excludeGeneIDs = []
							preGeneID = ""
							continue
						else:
							print "RF", "here is a stop"
							break
					currentSense = "-"
				scafPath.append(currentCtg)
				offset = gapStop
				preGeneID = nextGene.geneID
				currentVertex = nextVertex
				currentCtg = nameList[currentVertex]

				outstring = ">>>>(%d, %d): FF %d RR %d RF %d FR %d : %s \t%s" % (currentVertex, nextVertex, FF, RR, RF, FR, orientation, str(assemblyList))

			excludeGeneIDs = [preGeneID]
			scafPath.append(currentCtg)
			dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[currentCtg], dUpdateGFFs[scafName], scafName, offset, excludeGeneIDs)
			fOUTFASTA.write(">%s|%dbp|%s\n%s\n" % (scafName, len(sequence), ",".join(scafPath), sequence))
			fSCAFFPATH.write("%s\n" %(",".join(scafPath)))
			nCtgScaffolded += len(scafPath)
			print
	fSCAFFPATH.close()

	outScafGene = os.path.join(outDir, "%s.agouti.gene.contig_path" %(prefix))
	with open(outScafGene, 'w') as fOUTSCAFGENE:
		for k, v in sorted(dMergedGene2Ctgs.iteritems()):
			fOUTSCAFGENE.write(">%s\n%s\n" %(k, ','.join(v)))

	outGenePath = os.path.join(outDir, "%s.agouti.gene.gene_path" %(prefix))
	with open(outGenePath, 'w') as fOUTGENEPATH:
		for k, v in sorted(dMergedGene2Genes.iteritems()):
			fOUTGENEPATH.write(">%s\n%s\n" %(k, ','.join(v)))

	# other contigs need to be output
	sys.stderr.write("outputting contigs escaped from scaffolding ...\n")
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
	outGFF = os.path.join(outDir, "%s.agouti.gff" %(prefix))
	dFinalGFFs = dict(dGFFs, **dUpdateGFFs)
	numGenes = gff_outputter(dFinalGFFs, outGFF)

	sys.stderr.write("\n####Summary####\n")

	sys.stderr.write("number of contigs scaffoled: %d\n" %(nCtgScaffolded))
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
				fOUTGFF.write("%s\t%s\tExon\t%d\t%d\t.\t%s\t.\tID=%s.exon\n" %(k, geneModel.program, cdsStart, cdsStop, geneModel.strand, geneModel.geneID))
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
			print "loop", path
			loop2Remove.append(i)
			loop = 0
			for node in path:
				if node in visitedDict:
					del visitedDict[node]
	return loop2Remove, visitedDict

def get_gene_index(geneModels, curGeneID):
	print ">>>> get_gene_index", " ".join([x.geneID for x in geneModels]), "curGeneID", curGeneID
	currentGeneIndex = [x.geneID for x in geneModels].index(curGeneID)
	toLeft = currentGeneIndex
	toRight = len(geneModels) - 1 - currentGeneIndex
	return toLeft, toRight

def merge_gene_model(currentGene, nextGene, scafName,
					 numMergedGene, currentOffset, nextOffset):
	mergedGene = agGFF.AGOUTI_GFF()
	mergedGene.geneID = "AGOUTI_Merged_gene_%d" %(numMergedGene)
	mergedGene.geneStart = currentGene.geneStart + currentOffset
	mergedGene.geneStop = nextGene.geneStop + nextOffset
	mergedGene.program = "AGOUTI"
	mergedGene.ctgID = scafName
	if currentGene.strand == nextGene.strand:
		mergedGene.strand = currentGene.strand
	else:
		mergedGene.strand = '.'
	print ">>>> merge", "merged GeneID", mergedGene.geneID, "current LCDS", currentGene.lcds
	for i in range(len(currentGene.lcds)):
		currentGene.lcds[i] += currentOffset
	print ">>>> merge", "next LCDS before", nextGene.lcds
	for i in range(len(nextGene.lcds)):
		nextGene.lcds[i] += nextOffset
	print ">>>> merge", "next LCDS after", nextGene.lcds
	mergedGene.lcds = currentGene.lcds + nextGene.lcds
	print ">>>> merge", "finish merging", mergedGene.lcds
	return mergedGene

def update_gene_model(geneModels, updatedGeneModels,
					  newScafID, offset, excludes,
					  mergedGene=None):
	print ">>>> update"
	print ">>>> excludes to update", excludes
	tmpGeneModel = agGFF.AGOUTI_GFF()
	for i in range(len(geneModels)):
		tmpGeneModel = geneModels[i]
		if tmpGeneModel.geneID not in excludes:
			tmpGeneModel = geneModels[i]
			tmpGeneModel.ctgID = newScafID
			tmpGeneModel.geneStart += offset
			tmpGeneModel.geneStop += offset
			for j in range(len(tmpGeneModel.lcds)):
				tmpGeneModel.lcds[j] += offset
			updatedGeneModels.append(tmpGeneModel)
			print ">>>>", "update", tmpGeneModel.geneID, tmpGeneModel.lcds
	if mergedGene is not None:
		print ">>>> mergedGene to update", mergedGene.geneID
		updatedGeneModels.append(mergedGene)
		print ">>>>", "update", mergedGene.geneID, mergedGene.lcds
	updatedGeneIDs = []
	for i in range(len(updatedGeneModels)):
		updatedGeneIDs.append(updatedGeneModels[i].geneID)
	print "updated GeneIDs", updatedGeneIDs
	return updatedGeneModels, updatedGeneIDs

def reverse_gene_models(geneModels, lenCtg, excludeIDs=[]):
	'''
		reverse gene models
		call reverse_gene_model() function
	'''
	print ">>>>excludes to reverse", excludeIDs
	for i in xrange(len(geneModels)):
		if geneModels[i].geneID not in excludeIDs:
			geneModels[i] = reverse_gene_model(geneModels[i], lenCtg)
	return geneModels

def exclude_gene_model(geneModels, exclude):
	tmp = []
	for i in xrange(len(geneModels)):
		if geneModels[i].geneID != exclude:
			tmp.append(geneModels[i])
	return tmp

def reverse_gene_model(geneModel, lenCtg):
	'''
		reverse each gene, all CDS within the gene
	'''
	print ">>>>reverse, %d" %(lenCtg)
	tmpGeneStop =  lenCtg-geneModel.geneStart+1
	tmpGeneStart = lenCtg-geneModel.geneStop+1
	geneModel.geneStart = tmpGeneStart
	geneModel.geneStop = tmpGeneStop
	lcds = geneModel.lcds
	reverseLCDs = []
	print ">>>> before reverse", geneModel.geneID, "LCDS", lcds
	for j in range(len(lcds)-1, -1, -2):
		reverseLCDs += [lenCtg-lcds[j]+1, lenCtg-lcds[j-1]+1]
	print ">>>> after reverse", geneModel.geneID, "LCDS", reverseLCDs
	geneModel.lcds = reverseLCDs
	if geneModel.strand == '+':
		geneModel.strand = '-'
	else:
		geneModel.strand = '+'
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

def rc_seq(seq):
	"""
		return the reverse and complementary of
		the give sequence
	"""
	alphabetsInString = "ACGTNTGCANacgtntgcan"
	alphabets = { alphabetsInString[i]:alphabetsInString[i+5] for i
				  in range(20)
				  if i<5 or 10<=i<=14 }
	return "".join([alphabets[base] for base in seq[::-1]])
