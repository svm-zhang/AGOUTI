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

	# shut it off for now; working to improve it
	#agPathProgress.logger.info("[BEGIN] Checking consistency")
	#compare(dOriPaths, agoutiPaths, vertex2Name, outDir, prefix)
	#agPathProgress.logger.info("[DONE]")

	report_consistency(agoutiPaths, dOriPaths, vertex2Name, outDir, prefix)

	agPathProgress.logger.info("[BEGIN] Recovring original scaffolding")
	agoutiPaths, dCtgPair2GenePair, dSenses = recover_untouched_sequences(dOriPaths, agoutiPaths,
																		  vertex2Name, dGFFs,
																		  dCtgPair2GenePair,
																		  dSenses, agPathProgress,
																		  agPathDebug)
	agPathProgress.logger.info("[DONE]")

	return agoutiPaths, dCtgPair2GenePair, dSenses

def report_consistency(agoutiPaths, dOriPaths, vertex2Name, outDir, prefix, context=1):
	outConsist = os.path.join(outDir, "%s.consistency.nw" %(prefix))
	fOUT = open(outConsist, 'w')
	fOUT.write("source\tinteraction\ttarget\ttype\n")
	dPair = {}
	for agPath in agoutiPaths:
		agPath = [vertex2Name[k] for k in agPath]
		#print ">", agPath
		preCtg = agPath[0]
		for i in range(1, len(agPath)):
			curCtg = agPath[i]

			if len(preCtg.split('_')) == 1:
				preOrigScaf = preCtg
			else:
				preIndex = preCtg.split('_')[-1]
				preOrigScaf = preCtg.rstrip(preIndex).rstrip("_")
			if len(curCtg.split('_')) == 1:
				curOrigScaf = curCtg
			else:
				curIndex = curCtg.split('_')[-1]
				curOrigScaf = curCtg.rstrip(curIndex).rstrip("_")
			#print preCtg, preOrigScaf
			#print curCtg, curOrigScaf
			curIndex = int(curIndex)
			preIndex = int(preIndex)
			if preOrigScaf == curOrigScaf:
				if math.fabs(curIndex-preIndex) == 1:
					# consistent consecutive
					if (preCtg, curCtg) not in dPair and (curCtg, preCtg) not in dPair:
						fOUT.write("%s\t%s\t%s\tagouti_same\n"
								   %(preCtg, preOrigScaf, curCtg))
						dPair[preCtg, curCtg] = 1
						dPair[curCtg, preCtg] = 1
					#print "consistent consecutive"
				else:
					# consistent non-conseutive
					if (preCtg, curCtg) not in dPair and (curCtg, preCtg) not in dPair:
						fOUT.write("%s\t%s\t%s\tagouti_same\n"
								   %(preCtg, preOrigScaf, curCtg))
						dPair[preCtg, curCtg] = 1
						dPair[curCtg, preCtg] = 1
					tmpPreCtg = preCtg
					for j in range(preIndex+1, curIndex+1):
						tmpCurCtg = preOrigScaf + "_%d" %(j)
						if (tmpPreCtg, tmpCurCtg) not in dPair and \
						   (tmpCurCtg, tmpPreCtg) not in dPair:
							fOUT.write("%s\t%s\t%s\tskipped\n"
									   %(tmpPreCtg, preOrigScaf, tmpCurCtg))
							dPair[tmpPreCtg, tmpCurCtg] = 1
							dPair[tmpCurCtg, tmpPreCtg] = 1
						tmpPreCtg = tmpCurCtg
					#print "consistent non-consecutive"
			else:
				# inconsistent
				if (preCtg, curCtg) not in dPair and (curCtg, preCtg) not in dPair:
					fOUT.write("%s\t%s\t%s\tagouti_diff\n"
							   %(preCtg, preOrigScaf, curCtg))
					dPair[preCtg, curCtg] = 1
					dPair[curCtg, preCtg] = 1
				#print "inconsistent"
				start = preIndex-context
				stop = preIndex+context+1
				if start < 0:
					start = 0
				if stop > len(dOriPaths[preOrigScaf]):
					stop = len(dOriPaths[preOrigScaf])
				tmpPreCtg = preOrigScaf + "_%d" %(start)
				for j in range(start+1, stop):
					tmpCurCtg = preOrigScaf + "_%d" %(j)
					if (tmpPreCtg, tmpCurCtg) not in dPair and \
					   (tmpCurCtg, tmpPreCtg) not in dPair:
						fOUT.write("%s\t%s\t%s\toriginal\n"
								   %(tmpPreCtg, preOrigScaf, tmpCurCtg))
						dPair[tmpPreCtg, tmpCurCtg] = 1
						dPair[tmpCurCtg, tmpPreCtg] = 1
					tmpPreCtg = tmpCurCtg
				start = curIndex-context
				stop = curIndex+context+1
				if start < 0:
					start = 0
				if stop > len(dOriPaths[curOrigScaf]):
					stop = len(dOriPaths[curOrigScaf])
				tmpPreCtg = curOrigScaf + "_%d" %(start)
				for j in range(start+1, stop):
					tmpCurCtg = curOrigScaf + "_%d" %(j)
					if (tmpPreCtg, tmpCurCtg) not in dPair and \
					   (tmpCurCtg, tmpPreCtg) not in dPair:
						fOUT.write("%s\t%s\t%s\toriginal\n"
								   %(tmpPreCtg, curOrigScaf, tmpCurCtg))
						dPair[tmpPreCtg, tmpCurCtg] = 1
						dPair[tmpCurCtg, tmpPreCtg] = 1
					tmpPreCtg = tmpCurCtg
			preCtg = curCtg
	fOUT.close()
	#outAttr = os.path.join(outDir, "%s.consistency.attr" %(prefix))
	#fATTR = open(outAttr, 'w')
	#fATTR.write("node\ttype\n")
	#agNodes = {k:1 for path in agoutiPaths for k in path}
	#oriNodes = {v:1 for oriPath in dOriPaths.itervalues() for v in oriPath}
	#for oriNode in oriNodes.iterkeys():
	#	if oriNode not in agNodes:
	#		fATTR.write("%s\toriginal\n" %(oriNode))
	#	else:
	#		fATTR.write("%s\tagouti\n" %(oriNode))
	#fATTR.close()

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
		# len(pathLefts) == 0 means all contigs being touched by AGOUTI
		# len(pathLefts) == 1 means the contig is the original scaffold, no shred was done
		if len(pathLefts) <= 1:
			continue
		agPathDebug.debugger.debug("[RECOVER]\n")
		agPathDebug.debugger.debug("[RECOVER]\tpathLefts=%s" %(str(pathLefts)))
		agPathDebug.debugger.debug("[RECOVER]\tGetting consective path from pathLefts")
		recovPath = []
		preCtg = pathLefts[0]
		preCtgIndex = int(preCtg.split('_')[-1])
		tmpPath = [preCtg]
		for i in range(1, len(pathLefts)):
			curCtg = pathLefts[i]
			curCtgIndex = int(curCtg.split('_')[-1])
			agPathDebug.debugger.debug("[RECOVER]\tpreCtg=%s - preIndex=%d - curCtg=%s - curIndex=%d"
									   %(preCtg, preCtgIndex, curCtg, curCtgIndex))
			# this if/else makes sure only
			# consecutive contigs will be recovered
			if math.fabs(curCtgIndex-preCtgIndex) == 1:
				tmpPath.append(curCtg)
			else:
				if len(tmpPath) <= 1:
					tmpPath = [curCtg]
				else:
					recovPath.append(tmpPath)
					tmpPath = [curCtg]
			preCtg = curCtg
			preCtgIndex = int(preCtg.split('_')[-1])
		# add path with at least two contigs
		if len(tmpPath) > 1:
			recovPath.append(tmpPath)
		agPathDebug.debugger.debug("[RECOVER]\tGetting gene pair for each contig pair")
		for path in recovPath:
			agPathDebug.debugger.debug("[RECOVER]\t%s" %(str(path)))
			for i in range(1, len(path)):
				preCtg = path[i-1]
				curCtg = path[i]
				agPathDebug.debugger.debug("[RECOVER]\tpreCtg=%s - curCtg=%s"
										   %(preCtg, curCtg))
				preVertex = vertex2Name.index(preCtg)
				curVertex = vertex2Name.index(curCtg)
				preGeneID = ""
				curGeneID = ""
				preCtgGeneModel3 = None
				curCtgGeneModel5 = None
				# check if pre and cur contigs
				# have gene models on them
				if preCtg in dGFFs:
					preCtgGeneModel3 = dGFFs[preCtg][-1]
					if not preCtgGeneModel3.fake:
						tmpName = preCtgGeneModel3.geneID.split('_')
						preGeneID = tmpName[:-1]
				if curCtg in dGFFs:
					curCtgGeneModel5 = dGFFs[curCtg][0]
					if not curCtgGeneModel5.fake:
						tmpName = curCtgGeneModel5.geneID.split('_')
						curGeneID = tmpName[:-1]

				# the two contigs have no gene models
				# create fake0 genes and the merged gene won't be output
				if preGeneID == "" and curGeneID == "":
					# case 1
					preCtgGeneModel3 = create_fake0_gene(preCtg, m)
					dGFFs[preCtg].append(preCtgGeneModel3)
					m += 1
					agPathDebug.debugger.debug("[RECOVER]\t[CASE1]\t====preCtgGeneModels=%s"
											   %(str([k.geneID for k in dGFFs[preCtg]])))
					curCtgGeneModel5 = create_fake0_gene(curCtg, m)
					dGFFs[curCtg] = [curCtgGeneModel5] + dGFFs[curCtg]
					m += 1
					agPathDebug.debugger.debug("[RECOVER]\t[CASE1]\t====curCtgGeneModels=%s"
											   %(str([k.geneID for k in dGFFs[curCtg]])))
				# the two contigs have the same gene models
				# go ahead merging them
				# ignore if they are not consecutive
				elif preGeneID == curGeneID:
					# case 2
					n += 1
					agPathDebug.debugger.debug("[RECOVER]\t[CASE2]\t====preGeneID=%s - curGeneID=%s"
											   %(preGeneID, curGeneID))
				# the two contigs have two different gene models
				# The two genes cannot be merged
				# create fake0 gene model and the merged gene won't be output
				elif preGeneID != curGeneID:
					# case 3
					preCtgGeneModel3 = create_fake0_gene(preCtg, m)
					dGFFs[preCtg].append(preCtgGeneModel3)
					agPathDebug.debugger.debug("[RECOVER]\t[CASE3]\t====preCtgGeneModels=%s"
											   %(str([k.geneID for k in dGFFs[preCtg]])))
					m += 1
					curCtgGeneModel5 = create_fake0_gene(curCtg, m)
					dGFFs[curCtg] = [curCtgGeneModel5] + dGFFs[curCtg]
					m += 1
					agPathDebug.debugger.debug("[RECOVER]\t[CASE3]\t====curCtgGeneModels=%s"
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

def create_fake0_gene(ctg, m):
	geneModel = agGFF.AGOUTI_GFF()
	geneModel.setGene("%s_fakeO_%d" %(ctg, m),
							 0, 0, 1)
	return geneModel

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

