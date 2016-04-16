import sys
import os
import math

from lib import agouti_gff as agGFF
from lib import agouti_log as agLOG

def report_scaffold_path(paths, vertex2Name, outDir, prefix):
	'''
		output scaffolding path
		check conflicts with original paths
	'''
	outPathFile = os.path.join(outDir, "%s.agouti.scaffolding_paths.txt" %(prefix))
	with open(outPathFile, 'w') as fSCAFPATH:
		for i, path in enumerate(paths):
			scafName = prefix + "_scaf_%d" %(i+1)
			fSCAFPATH.write(">%s\n%s\n" %(scafName, ",".join([vertex2Name[k] for k in path])))

def agouti_path_main(agoutiPaths, dSenses, vertex2Name,
					 dGFFs, dCtgPair2GenePair,
					 oriScafPathFile, outDir, prefix):
	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_path")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)
	agPathProgress = agLOG.PROGRESS_METER(moduleName)
	agPathProgress.logger.info("Analyzing scaffolding paths")
	outDebugFile = os.path.join(moduleOutDir, prefix) + ".agouti_path.debug"
	agPathDebug = agLOG.DEBUG("SHREDDER", outDebugFile)
	agPathProgress.logger.info("[BEGIN] Reading file with shred info")
	dOriPaths, dOriGaps = read_original_path(oriScafPathFile, agPathProgress)
	agPathProgress.logger.info("[DONE]")
	agPathProgress.logger.info("[BEGIN] Checking consistency")
	compare(dOriPaths, agoutiPaths, vertex2Name, outDir, prefix)
	agPathProgress.logger.info("[DONE]")
	agPathProgress.logger.info("[BEGIN] Recovring original scaffolding")
	agoutiPaths, dCtgPair2GenePair, dSenses = recover_untouched_sequences(dOriPaths, agoutiPaths,
																 vertex2Name, dGFFs,
																 dCtgPair2GenePair,
																 dSenses, agPathProgress,
																 agPathDebug)
	agPathProgress.logger.info("[DONE]")

	return agoutiPaths, dCtgPair2GenePair, dSenses

def recover_untouched_sequences(dOriPaths, agoutiPaths, vertex2Name,
								dGFFs, dCtgPair2GenePair, dSenses,
								agPathProgress, agPathDebug):
	'''
		Recover scaffolds from contigs untouched by AGOUTI
	'''

	agPathDebug.debugger.debug("[RECOVER]\t####")
	agPathDebug.debugger.debug("[RECOVER]\tNumber of paths AGOUTI scaffolded: %d"
							   %(sum([len(path) for path in agoutiPaths])))
	agPathDebug.debugger.debug("[RECOVER]\tNumber of gene models in dGFFs: %d"
							   %(sum([len(geneModels) for geneModels in dGFFs.itervalues()])))
	agPathDebug.debugger.debug("[RECOVER]\tNumber of fake gene models in dGFFs: %d"
							   %(len([geneModel
									  for geneModels in dGFFs.itervalues()
									  for geneModel in geneModels
									  if  geneModel.fake == 1])))
	agPathDebug.debugger.debug("[RECOVER]\t####")
	# get a dictionary of contigs being scaffolded
	scaffoldedCtgs = {vertex2Name[k]:0 for path in agoutiPaths for k in path}
	n = 0				# number of canonical joining gene models
	m = 0				# number of fakeO gene model created
	tmpPaths = []
	for k, path in dOriPaths.iteritems():
		pathLefts = [ctg for ctg in path if ctg not in scaffoldedCtgs]
		if not pathLefts:
			continue
		agPathDebug.debugger.debug("[RECOVER]\tpathLefts=%s" %(str(pathLefts)))
		agPathDebug.debugger.debug("[RECOVER]\tGetting consective path from pathLefts")
		recovPath = []
		preCtg = pathLefts[0]
		preCtgIndex = int(preCtg.split('_')[1])
		tmpPath = [preCtg]
		for i in range(1, len(pathLefts)):
			curCtg = pathLefts[i]
			curCtgIndex = int(curCtg.split('_')[1])
			agPathDebug.debugger.debug("[RECOVER]\tpreCtg=%s - preIndex=%d - curCtg=%s - curIndex=%d"
									   %(preCtg, preCtgIndex, curCtg, curCtgIndex))
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
		agPathDebug.debugger.debug("[RECOVER]\tGetting gene pair for each contig pair")
		for path in recovPath:
			for i in range(1, len(path)):
				preCtg = path[i-1]
				curCtg = path[i]
				agPathDebug.debugger.debug("[RECOVER]\tpreCtg=%s - curCtg=%s"
										   %(preCtg, curCtg))
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

				agPathDebug.debugger.debug("[RECOVER]\t====preGeneIndex=%d - curGeneIndex=%d"
										   %(preIndex, curIndex))
				# both index should not be -1 if they are the
				# joining gene model between the two contigs
				if preIndex != -1 and curIndex != -1:
					# and they have to be consecutive
					if math.fabs(curIndex - preIndex) == 1:
						n += 1
						agPathDebug.debugger.debug("[RECOVER]\t====preGeneID=%s - curGeneID=%s"
												   %(preGeneID, curGeneID))
					# create fakeO genes as joining gene models
					else:
						preCtgGeneModel3 = agGFF.AGOUTI_GFF()
						preCtgGeneModel3.setGene("%s_fakeO_%d" %(preCtg, m),
												 0, 0, 1)
						dGFFs[preCtg].append(preCtgGeneModel3)
						m += 1
						agPathDebug.debugger.debug("[RECOVER\t====preCtgGeneModels=%s"
												   %(str([k.geneID for k in dGFFs[preCtg]])))
						curCtgGeneModel5 = agGFF.AGOUTI_GFF()
						curCtgGeneModel5.setGene("%s_fakeO_%d" %(curCtg, m),
												 0, 0, 1)
						dGFFs[curCtg] = [curCtgGeneModel5] + dGFFs[curCtg]
						m += 1
						agPathDebug.debugger.debug("[RECOVER\t====curCtgGeneModels=%s"
												   %(str([k.geneID for k in dGFFs[curCtg]])))
				# create fakeO gnes as joining gene models
				else:
					preCtgGeneModel3 = agGFF.AGOUTI_GFF()
					preCtgGeneModel3.setGene("%s_fakeO_%d" %(preCtg, m),
											 0, 0, 1)
					dGFFs[preCtg].append(preCtgGeneModel3)
					m += 1
					agPathDebug.debugger.debug("[RECOVER\t====preCtgGeneModels=%s"
											   %(str([k.geneID for k in dGFFs[preCtg]])))
					curCtgGeneModel5 = agGFF.AGOUTI_GFF()
					curCtgGeneModel5.setGene("%s_fakeO_%d" %(curCtg, m),
											 0, 0, 1)
					dGFFs[curCtg] = [curCtgGeneModel5] + dGFFs[curCtg]
					m += 1
					agPathDebug.debugger.debug("[RECOVER\t====curCtgGeneModels=%s"
											   %(str([k.geneID for k in dGFFs[curCtg]])))
				# record joining gene models for the two contigs
				dCtgPair2GenePair[preVertex, curVertex] = [preCtgGeneModel3, curCtgGeneModel5]
				agPathDebug.debugger.debug("[RECOVER]\t====preID=%s preCDS=%s"
										   %(preCtgGeneModel3.geneID, str(preCtgGeneModel3.lcds)))
				agPathDebug.debugger.debug("[RECOVER]\t====curID=%s curCDS=%s"
										   %(curCtgGeneModel5.geneID, str(curCtgGeneModel5.lcds)))
				dSenses[preVertex, curVertex] = [('+', '-')]
		for path in recovPath:
			# convert from name to vertex
			agoutiPaths.append([vertex2Name.index(k) for k in path])
	agPathDebug.debugger.debug("[RECOVER]\t####")
	agPathDebug.debugger.debug("[RECOVER]\tNumber of paths to update: %d"
							   %(sum([len(path) for path in agoutiPaths])))
	agPathDebug.debugger.debug("[RECOVER]\tNumber of gene models in dGFFs: %d"
							   %(sum([len(geneModels) for geneModels in dGFFs.itervalues()])))
	agPathDebug.debugger.debug("[RECOVER]\tNumber of fake gene models in dGFFs: %d"
							   %(len([geneModel
									  for geneModels in dGFFs.itervalues()
									  for geneModel in geneModels
									  if  geneModel.fake == 1])))
	agPathDebug.debugger.debug("[RECOVER]\tNumber of fakeO gene models in dGFFs: %d"
							   %(len([geneModel
									  for geneModels in dGFFs.itervalues()
									  for geneModel in geneModels
									  if  geneModel.fake == 1
										  and geneModel.geneStop == 0])))
	agPathDebug.debugger.debug("[RECOVER]\t####")
	agPathProgress.logger.info(("Number of pairs of shredded contigs "
								"having joining gene models: %d" %(n)))
	return agoutiPaths, dCtgPair2GenePair, dSenses

def read_original_path(oriScafPathFile, agPathProgress):
	'''
		Read original scaffolding path
	'''
	dOriPaths = {}
	dOriGaps = {}
	oriScafID = ""
	with open(oriScafPathFile, 'r') as fIN:
		oriPath = []
		nLines = 0
		for line in fIN:
			nLines += 1
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
						agPathProgress.logger.error("Error: ctg1=%s ctg2=%s occurred twice")
						agPathProgress.logger.error("Please check file: %s at line %d"
													%(oriScafPathFile, nLines))
						agPathProgress.logger.error(("The error probably comes from incorrect"
													 "recording shred info"))
						agPathProgress.logger.error("Please report this bug")
						sys.exit(1)
					if len(oriPath) == 0:
						oriPath += [curCtg, nexCtg]
					else:
						oriPath.append(nexCtg)
				else:
					oriPath.append(curCtg)
		if oriPath:
			dOriPaths[oriScafID] = oriPath
	agPathProgress.logger.info("Number of shredded contigs parsed: %d"
								%(sum([len(v) for v in dOriPaths.itervalues()])))

	return dOriPaths, dOriGaps

def compare(dOriPaths, agPaths, vertex2Name, outDir, prefix):
	fCONFLICT = None
	if dOriPaths:
		outConflicts = os.path.join(outDir, "%s.agouti_vs_original.compare.txt" %(prefix))
		fCONFLICT = open(outConflicts, 'w')
		if dOriPaths:
			for i, path in enumerate(agPaths):
				path = [vertex2Name[k] for k in path]
				conflictType = check_consistency(dOriPaths, path)
				scafName = "%s_scaf_%d" %(prefix, i)
				if conflictType is not None:
					fCONFLICT.write("%s\t%s\t%s\n" %(scafName, conflictType,
													 ",".join(path)))

def check_consistency(dOriPaths, agPath):
	for i in range(1, len(agPath)):
		preCtg = agPath[i-1]
		curCtg = agPath[i]
		preCtgScafID = preCtg.split('_')[0]
		curCtgScafID = curCtg.split('_')[0]
		# cases where no shred on the original scaffold
		if preCtgScafID == preCtg or curCtgScafID == curCtg:
			return "INTERSCAFFOLDING"
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

