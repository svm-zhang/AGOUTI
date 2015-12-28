import sys
import os
import math

from lib import agouti_gff as agGFF

def report_scaffold_path(paths, vertex2Name, outDir, prefix):
	'''
		output scaffolding path
		check conflicts with original paths
	'''
	outPathFile = os.path.join(outDir, "%s.agouti.scaffolding_paths.txt" %(prefix))
	with open(outPathFile, 'w') as fSCAFPATH:
		for i, path in enumerate(paths):
			scafName = prefix + "_scaf_%d" %(i)
			fSCAFPATH.write(">%s\n%s\n" %(scafName, ",".join([vertex2Name[k] for k in path])))

def agouti_path_main(agoutiPaths, dSenses, vertex2Name,
					 dGFFs, dCtgPair2GenePair,
					 oriScafPathFile, outDir, prefix):
	print "Analyzing scaffolding paths"
	dOriPaths, dOriGaps = read_original_path(oriScafPathFile)
	agoutiPaths, dCtgPair2GenePair, dSenses = recover_untouched_sequences(dOriPaths, agoutiPaths,
																 vertex2Name, dGFFs,
																 dCtgPair2GenePair,
																 dSenses)

	return agoutiPaths, dCtgPair2GenePair, dSenses

def recover_untouched_sequences(dOriPaths, agoutiPaths, vertex2Name,
								dGFFs, dCtgPair2GenePair, dSenses):
	# get a dictionary of contigs being scaffolded
	scaffoldedCtgs = {vertex2Name[k]:0 for path in agoutiPaths for k in path}

	#print len(agoutiPaths)
	#print sum([len(path) for path in agoutiPaths])
	#print sum([len(geneModels) for k, geneModels in dGFFs.iteritems()])
	#print len([geneModel for k, geneModels in dGFFs.iteritems() for geneModel in geneModels if geneModel.fake == 0])
	#print len([geneModel for k, geneModels in dGFFs.iteritems() for geneModel in geneModels if geneModel.fake == 1])
	n = 0	# number of canonical joining gene models
	m = 0	# number of fakeO gene model created
	tmpPaths = []
	for k, path in dOriPaths.iteritems():
		pathLefts = [ctg for ctg in path if ctg not in scaffoldedCtgs]
		if not pathLefts:
			continue
		#print "pathLefts", pathLefts
		recovPath = []
		preCtg = pathLefts[0]
		preCtgIndex = int(preCtg.split('_')[1])
		tmpPath = [preCtg]
		for i in range(1, len(pathLefts)):
			curCtg = pathLefts[i]
			curCtgIndex = int(curCtg.split('_')[1])
			#print preCtg, preCtgIndex, curCtg, curCtgIndex
			if math.fabs(curCtgIndex-preCtgIndex) == 1:
				tmpPath.append(curCtg)
			else:
				# add path with at least two contigs
				if len(tmpPath) <= 1:
					tmpPath = [curCtg]
				else:
					recovPath.append(tmpPath)
					tmpPath = [curCtg]
			preCtg = curCtg
			preCtgIndex = int(preCtg.split('_')[1])
		# add path with at least two contigs
		if len(tmpPath) > 1:
			recovPath.append(tmpPath)
		for path in recovPath:
			for i in range(1, len(path)):
				preCtg = path[i-1]
				curCtg = path[i]
				print "preCtg", preCtg, "curCtg", curCtg
				preVertex = vertex2Name.index(preCtg)
				curVertex = vertex2Name.index(curCtg)
				preIndex = -1
				curIndex = -1
				preCtgGeneModel3 = None
				curCtgGeneModel5 = None
				# check if pre and cur contigs
				# have gene models on them
				if preCtg in dGFFs:
					preCtgGeneModel3 = dGFFs[preCtg][-1]
					if not preCtgGeneModel3.fake:
						preGeneID = preCtgGeneModel3.geneID
						preIndex = int(preGeneID.split('_')[1])
				if curCtg in dGFFs:
					curCtgGeneModel5 = dGFFs[curCtg][0]
					if not curCtgGeneModel5.fake:
						curGeneID = curCtgGeneModel5.geneID
						curIndex = int(curGeneID.split('_')[1])

				print "preIndex", preIndex, "curIndex", curIndex
				# both index should not be -1 if they are the
				# joining gene model between the two contigs
				if preIndex != -1 and curIndex != -1:
					# and they have to be consecutive
					if math.fabs(curIndex - preIndex) == 1:
						n += 1
						print preGeneID, curGeneID
					# create fakeO genes as joining gene models
					else:
						preCtgGeneModel3 = agGFF.AGOUTI_GFF()
						preCtgGeneModel3.setGene("%s_fakeO_%d" %(preCtg, m),
												 0, 0, 1)
						dGFFs[preCtg].append(preCtgGeneModel3)
						m += 1
						curCtgGeneModel5 = agGFF.AGOUTI_GFF()
						curCtgGeneModel5.setGene("%s_fakeO_%d" %(curCtg, m),
												 0, 0, 1)
						dGFFs[curCtg] = [curCtgGeneModel5] + dGFFs[curCtg]
						m += 1
				# create fakeO gnes as joining gene models
				else:
					preCtgGeneModel3 = agGFF.AGOUTI_GFF()
					preCtgGeneModel3.setGene("%s_fakeO_%d" %(preCtg, m),
											 0, 0, 1)
					dGFFs[preCtg].append(preCtgGeneModel3)
					m += 1
					print [k.geneID for k in dGFFs[preCtg]]
					curCtgGeneModel5 = agGFF.AGOUTI_GFF()
					curCtgGeneModel5.setGene("%s_fakeO_%d" %(curCtg, m),
											 0, 0, 1)
					dGFFs[curCtg] = [curCtgGeneModel5] + dGFFs[curCtg]
					m += 1
					print [k.geneID for k in dGFFs[curCtg]]
				# record joining gene models for the two contigs
				dCtgPair2GenePair[preVertex, curVertex] = [preCtgGeneModel3, curCtgGeneModel5]
				print "preID", preCtgGeneModel3.geneID, preCtgGeneModel3.lcds
				print "curID", curCtgGeneModel5.geneID, curCtgGeneModel5.lcds
				dSenses[preVertex, curVertex] = [('+', '-')]
		for path in recovPath:
			# convert from name to vertex
			agoutiPaths.append([vertex2Name.index(k) for k in path])
	#print sum([len(path) for path in tmpPaths])
	#print len(tmpPaths)
	#print sum([len(geneModels) for k, geneModels in dGFFs.iteritems()])
	#print len([geneModel for k, geneModels in dGFFs.iteritems() for geneModel in geneModels if geneModel.fake == 0])
	#print len([geneModel for k, geneModels in dGFFs.iteritems() for geneModel in geneModels if geneModel.fake == 1])
	print "Number of pairs of shredded contigs having joining gene models: %d" %(n)
	return agoutiPaths, dCtgPair2GenePair, dSenses

def read_original_path(oriScafPathFile):
	dOriPaths = {}
	dOriGaps = {}
	oriScafID = ""
	with open(oriScafPathFile, 'r') as fIN:
		oriPath = []
		for line in fIN:
			if line.startswith('>'):
				if len(oriPath) > 0:
					dOriPaths[oriScafID] = oriPath
					oriPath = []
				oriScafID = line.strip()[1:]
			else:
				tmpLine = line.strip().split("\t")
				curCtg = tmpLine[0]
				nexCtg = tmpLine[1]
				if tmpLine[1] != "NA":
					gapLen = int(tmpLine[2])
					if (curCtg, nexCtg) not in dOriGaps:
						dOriGaps[curCtg, nexCtg] = gapLen
					else:
						print "Error: %s %s occurred already"
						print "Check the give shred info file"
						sys.exit(1)
					if len(oriPath) == 0:
						oriPath += [curCtg, nexCtg]
					else:
						oriPath.append(nexCtg)
				else:
					oriPath.append(curCtg)
		if oriPath:
			dOriPaths[oriScafID] = oriPath

	return dOriPaths, dOriGaps

def compare():
#	fCONFLICT = None
#	if dOriPaths:
#		outConflicts = os.path.join(outDir, "%s.agouti_vs_original.compare.txt" %(prefix))
#		fCONFLICT = open(outConflicts, 'w')
#			if dOriPaths:
#				conflictType = agSUM.check_consistency(dOriPaths, path)
#				if conflictType is not None:
#					fCONFLICT.write("%s\t%s\t%s\n" %(scafName, conflictType,
#													 ",".join(path)))
	pass

def check_consistency(dOriPaths, agPath):
	for i in range(1, len(agPath)):
		preCtg = agPath[i-1]
		curCtg = agPath[i]
		preCtgScafID = preCtg.split('_')[0]
		curCtgScafID = curCtg.split('_')[0]
		if preCtgScafID != curCtgScafID:
			preCtgScafIndex = dOriPaths[preCtgScafID].index(preCtg)
			curCtgScafIndex = dOriPaths[curCtgScafID].index(curCtg)
			preNCtg = len(dOriPaths[preCtgScafID])
			curNCtg = len(dOriPaths[curCtgScafID])
			if ((preCtgScafIndex == 0 or preCtgScafIndex == preNCtg-1) and
				(curCtgScafIndex == 0 or curCtgScafIndex == curNCtg-1)):
				return "INTERSCAFFOLDING"
			else:
				return "CONFLICT"
	agStartCtg = agPath[0]
	oriScafID = agStartCtg.split('_')[0]
	if oriScafID in dOriPaths:
		oriPath = dOriPaths[oriScafID]
		order = -1
		for i in range(1, len(agPath)):
			agCtgPre = agPath[i-1]
			agCtgScafIndexPre = int(agCtgPre.split('_')[1])
			agCtgNext = agPath[i]
			agCtgScafIndexNext = int(agCtgNext.split('_')[1])
			jump = agCtgScafIndexNext - agCtgScafIndexPre
			if order == -1:
				if jump <= -2:
					return "NONCONSECUTIVE"
				elif jump == -1:
					order = 1
				elif jump >= 2:
					return "NONCONSECUTIVE"
				elif jump == 1:
					order = 0
			else:
				if math.fabs(jump) >= 2:
					return "NONCONSECUTIVE"
				else:
					if order == 0 and jump < 0:
						return "SWITCHORDER"
					elif order == 1 and jump > 0:
						return "SWITCHORDER"
	return

