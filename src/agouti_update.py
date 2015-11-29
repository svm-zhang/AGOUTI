import os
import sys
import collections

from lib import agouti_log as agLOG
from lib import agouti_gff as agGFF
from src import agouti_sequence as agSeq
from src import agouti_summarize as agSUM

def agouti_update(pathList, dSeqs, seqNames,
				  edgeSenseDict, dGFFs,
				  dCtgPair2GenePair, outDir, prefix,
				  oriScafPathFile,
				  nFills=1000, debug=0):

	moduleName = os.path.basename(__file__).split('.')[0].upper()
	moduleOutDir = os.path.join(outDir, "agouti_update")
	if not os.path.exists(moduleOutDir):
		os.makedirs(moduleOutDir)

	progressLogFile = os.path.join(moduleOutDir, "%s.agouti_update.progressMeter" %(prefix))
	global agUPDATEProgress
	agUPDATEProgress = agLOG.PROGRESS_METER(moduleName)
	agUPDATEProgress.add_file_handler(progressLogFile)
	if debug:
		debugLogFile = os.path.join(moduleOutDir, "%s.agouti_update.debug" %(prefix))
		global agUPDATEDebug
		agUPDATEDebug = agLOG.DEBUG(moduleName, debugLogFile)

	agUPDATEProgress.logger.info("[BEGIN] Updating gene models")
	scafID = 0

	outFasta = os.path.join(outDir, "%s.agouti.fasta" %(prefix))
	fFASTA = open(outFasta, 'w')
	dUpdateGFFs = collections.defaultdict(list)
	dMergedGene2Ctgs = collections.defaultdict(list)
	dMergedGene2Genes = collections.defaultdict(list)
	numMergedGene = 0
	nCtgScaffolded = 0
	scaffoldedCtgs = {}
	seqLens = []
	agPaths = []
	mergedGenes = []
	for i in range(len(pathList)):
		path = pathList[i]
		scafID += 1
		scafName = prefix + "_scaf_%d" %(scafID)

		curVertex = path[0]
		sequence = dSeqs[curVertex]
		curSense = "+"
		curCtg = seqNames[curVertex]
		preCtg = ""
		scafPath = []
		assemblyList = curVertex
		preGeneID, curGeneID = "", ""
		mergedGene = agGFF.AGOUTI_GFF()
		preMergedGene = agGFF.AGOUTI_GFF()
		gapStart, gapStop = 0, 0
		offset = 0
		orientation = ""
		updatedGeneIDs = []
		mergedGenesPerPath = []
		excludeGeneIDs = []
		for nextVertex in path[1:]:
			nextCtg = seqNames[nextVertex]

			if preCtg == "":
				if debug:
					agUPDATEDebug.debugger.debug("UPDATE_MAIN\t>scaf_%d - path - %s"
												 %(scafID,
												  str([seqNames[vertex] for vertex in path])))
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tcurVertex - %d - %s - nextVertex - %d - %s"
											 %(curVertex, curCtg, nextVertex, nextCtg))

			#currentGene, nextGene = ctgpair2genepair(dCtgPair2GenePair, curCtg, nextCtg)
			currentGene, nextGene = ctgpair2genepair(dCtgPair2GenePair, curVertex, nextVertex)
			#!!! I should not break here, should continue#
			if currentGene is None and nextGene is None:
				agUPDATEProgress.error("%s - %s found no gene models joining them"
									   %(curCtg, nextCtg))
				agUPDATEProgress.error("This is NOT EXPECTED, REPORT!")
				sys.exit(1)
			curGeneID = currentGene.geneID
			excludeGeneIDs = [preGeneID] + [curGeneID]
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafName - %s" %(scafName))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tpreGene - %s - curGene - %s - nextGene - %s"
											 %(preGeneID, currentGene.geneID, nextGene.geneID))

			FR, FF, RR, RF = get_orientation_counts(curVertex, nextVertex, edgeSenseDict)
			if curSense == "-":
				temp1 = FR
				temp2 = FF
				FR = RR
				FF = RF
				RR = temp1
				RF = temp2
			orientation = decide_orientation(FR, FF, RR, RF)

			gapStart = gapStop + 1 + len(dSeqs[curVertex])
			gapStop = gapStart + nFills - 1
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tcurSense - %s" %(curSense))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tFF - %d - RR - %d - RF - %d - FR - %d"
											 %(FF, RR, RF, FR))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\toffset - %d - curCtgLen - %d"
											 %(offset, len(dSeqs[curVertex])))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tgapstart - %d - gapstop - %d"
											 %(gapStart, gapStop+1))
			valid = 0
			if orientation == "FR":
				if curGeneID !=  preGeneID:
					numMergedGene += 1
					mergedGene = merge_gene_model(currentGene, nextGene, scafName,
												  numMergedGene, offset, gapStop,
												  debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
					dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																			  scafName, offset, excludeGeneIDs,
																			  debug, mergedGene)
				else:
					mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
												  numMergedGene, 0, gapStop, debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
					indexMerged = updatedGeneIDs.index(mergedGene.geneID)
					dUpdateGFFs[scafName][indexMerged] = mergedGene
#				if debug:
#					agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tmergedGene - %s - dMergedGene2Ctgs - %s"
#												 %(mergedGene.geneID, str(dMergedGene2Ctgs[mergedGene.geneID])))
				preMergedGene = mergedGene
				assemblyList = ((assemblyList, "+"), (nextVertex, "+"))
				sequence += 'N'*nFills + dSeqs[nextVertex]
				curSense = "+"
			elif orientation == "FF":
				nextGene = reverse_gene_model(nextGene, len(dSeqs[nextVertex]), debug)
				if curGeneID !=  preGeneID:
					numMergedGene += 1
					mergedGene = merge_gene_model(currentGene, nextGene, scafName,
												  numMergedGene, offset, gapStop,
												  debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
					dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																			  scafName, offset, excludeGeneIDs,
																			  debug, mergedGene)
				else:
					mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
												  numMergedGene, 0, gapStop, debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
					indexMerged = updatedGeneIDs.index(mergedGene.geneID)
					dUpdateGFFs[scafName][indexMerged] = mergedGene
				preMergedGene = mergedGene
				assemblyList = ((assemblyList, "+"), (nextVertex, "-"))
				sequence += 'N'*nFills + agSeq.rc_seq(dSeqs[nextVertex])
				curSense = "-"
			elif orientation == "RR":
				if currentGene.geneID != preGeneID:
					dGFFs[curCtg] = reverse_gene_models(dGFFs[curCtg], len(dSeqs[curVertex]), debug)
					numMergedGene += 1
					mergedGene = merge_gene_model(currentGene, nextGene, scafName,
												  numMergedGene, offset, gapStop,
												  debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
					dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																			  scafName, offset, excludeGeneIDs,
																			  debug, mergedGene)
				else:
					dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																			  scafName, offset, excludeGeneIDs,
																			  debug)
					dUpdateGFFs[scafName] = reverse_gene_models(dUpdateGFFs[scafName], gapStart-1, debug)
					mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
												  numMergedGene, 0, gapStop,
												  debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
					indexMerged = updatedGeneIDs.index(mergedGene.geneID)
					dUpdateGFFs[scafName][indexMerged] = mergedGene
				assemblyList = ((assemblyList, "-"), (nextVertex, "+"))
				sequence = agSeq.rc_seq(sequence) + \
						   'N'*nFills + dSeqs[nextVertex]
				preMergedGene = mergedGene
					#print ">>>> RR zone preMergedGene", preMergedGene.geneID, preMergedGene.geneStart, preMergedGene.geneStop
				curSense = "+"
			elif orientation == "RF":
				if currentGene.geneID != preGeneID:
					dGFFs[curCtg] = reverse_gene_models(dGFFs[curCtg], len(dSeqs[curVertex]), debug)
					nextGene = reverse_gene_model(nextGene, len(dSeqs[nextVertex]), debug)
					numMergedGene += 1
					mergedGene = merge_gene_model(currentGene, nextGene, scafName,
												  numMergedGene, offset, gapStop,
												  debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
					dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																			  scafName, offset, excludeGeneIDs,
																			  debug, mergedGene)
				else:
					dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																			  scafName, offset, excludeGeneIDs,
																			  debug)
					dUpdateGFFs[scafName] = reverse_gene_models(dUpdateGFFs[scafName],
																gapStop+len(dSeqs[curVertex]),
																debug)
					mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
												  numMergedGene, 0, gapStop,
												  debug)
					dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
					dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
					indexMerged = updatedGeneIDs.index(mergedGene.geneID)
					dUpdateGFFs[scafName][indexMerged] = mergedGene
				assemblyList = ((assemblyList, "-"), (nextVertex, "-"))
				sequence = agSeq.rc_seq(sequence) + \
						   'N'*nFills + \
						   agSeq.rc_seq(dSeqs[nextVertex])
				preMergedGene = mergedGene
					#print ">>>> RF zone preMergedGene", preMergedGene.geneID, preMergedGene.geneStart, preMergedGene.geneStop
				curSense = "-"
			scafPath.append(curCtg)
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafPath updates- %s"
											 %(str(scafPath)))
			mergedGenesPerPath.append(mergedGene.geneID)
			offset = gapStop
			preGeneID = nextGene.geneID
			preCtg = curCtg
			curVertex = nextVertex
			curCtg = seqNames[curVertex]

		excludeGeneIDs = [preGeneID]
		if debug:
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tappend last curCtg - %s" %(curCtg))
		scafPath.append(curCtg)
		if debug:
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafPath - %s"
										 %(str(scafPath)))
		mergedGenes.append(mergedGenesPerPath)
		dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg],
																  dUpdateGFFs[scafName],
																  scafName, offset,
																  excludeGeneIDs, debug)
		fFASTA.write(">%s|%dbp|%s\n%s\n"
						%(scafName, len(sequence), ",".join(scafPath), sequence))
		seqLens.append(len(sequence))
		agPaths.append(scafPath)
		nCtgScaffolded += len(scafPath)
		scaffoldedCtgs.update(dict((contig, 1) for contig in scafPath))
		if debug:
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tupdatedGeneIDs - %s"
										 %(str(updatedGeneIDs)))
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tmergedGenesPerPath - %s"
										 %(str(mergedGenesPerPath)))
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t-------------------------------------")

	dOriPaths = None
	dOriGaps = None
	if oriScafPathFile:
		dOriPaths, dOriGaps = agSUM.read_original_path(oriScafPathFile)
	# other contigs need to be output
	agUPDATEProgress.logger.info("Finalizing sequences")
	numLeft = 0
	for vertex in dSeqs:
		if seqNames[vertex] in scaffoldedCtgs:
			if seqNames[vertex] in dGFFs:
				del dGFFs[seqNames[vertex]]
			continue
		numLeft += 1
		fFASTA.write(">%s\n%s\n" % (seqNames[vertex], dSeqs[vertex]))
		seqLens.append(len(dSeqs[vertex]))
	fFASTA.close()
	n50 = agSeq.get_assembly_NXX(seqLens)

	agUPDATEProgress.logger.info("Summarizing AGOUTI scaffolding paths")
	summarize_scaffold_path(agPaths, outDir, prefix, dOriPaths)

	agUPDATEProgress.logger.info("Summarizing AGOUTI gene paths")
	summarize_gene_path(dMergedGene2Genes, dMergedGene2Ctgs,
						outDir, prefix)

	agUPDATEProgress.logger.info("Outputting updated Gene Moddels")
	for vertex in dSeqs:
		if seqNames[vertex] in scaffoldedCtgs:
			if seqNames[vertex] in dGFFs:
				del dGFFs[seqNames[vertex]]
	dFinalGFFs = dict(dGFFs, **dUpdateGFFs)
	numGenes = output_gff(dFinalGFFs, dMergedGene2Ctgs, dMergedGene2Genes,
						  outDir, prefix)

	agUPDATEProgress.logger.info("-----------Summary-----------")

	agUPDATEProgress.logger.info("number of contigs scaffoled: %d" %(nCtgScaffolded))
	agUPDATEProgress.logger.info("number of scaffolds: %d" %(scafID))
#	agUPDATEProgress.logger.info("number of contigs found no links: %d" %(numLeft))
	agUPDATEProgress.logger.info("number of contigs in the final assembly: %d" %(len(seqLens)))
	agUPDATEProgress.logger.info("Final assembly N50: %d" %(n50))
	agUPDATEProgress.logger.info("Final number of genes: %d" %(numGenes))
	agUPDATEProgress.logger.info("Succeeded")

def recover_original_scaffold():
	'''
		not finish
	'''
	if dOriPaths is None:
		for vertex in sorted(dSeqs):
			if seqNames[vertex] not in scaffoldedCtgs:
				fFASTA.write(">%s\n%s\n" % (seqNames[vertex], dSeqs[vertex]))
				seqLens.append(len(dSeqs[vertex]))
	else:
		for _, oriPath in dOriPaths.iteritems():
			if len(oriPath) == 1:
				index = seqNames.index(oriPath[0])
				fFASTA.write(">%s\n%s\n" % (seqNames[index], dSeqs[index]))
				seqLens.append(len(dSeqs[index]))
				continue
			print oriPath
			untouchedCtgs = [k for k in oriPath if k not in scaffoldedCtgs]
			print untouchedCtgs
			sequence = ""
			preCtg = untouchedCtgs[0]
			preIndex = seqNames.index(preCtg)
			sequence = dSeqs[preIndex]
			for i in range(1, len(untouchedCtgs)):
				curCtg = untouchedCtgs[i]
				if (preCtg, curCtg) in dOriGaps:
					gapLen = dOriGaps[preCtg, curCtg]
				else:
					fFASTA.write(">%s\n%s\n" % ("haha_need_to_fix", sequence))
					seqLens.append(len(sequence))
					sequence = ""
					preCtg = curCtg
					continue
				preIndex = seqNames.index(preCtg)
				curIndex = seqNames.index(curCtg)
				sequence += gapLen * 'N' + \
							dSeqs[curIndex]
				print "preCtg", preCtg, "curCtg", curCtg
				print len(sequence)
				preCtg = curCtg
			if sequence:
				fFASTA.write(">%s\n%s\n" % ("haha_need_to_fix", sequence))
				seqLens.append(len(sequence))
			sys.exit()

def summarize_scaffold_path(agPaths, outDir, prefix, dOriPaths):
	'''
		output scaffolding path
		check conflicts with original paths
	'''
	outScaffPath = os.path.join(outDir, "%s.agouti.scaffolding_paths.txt" %(prefix))
	fCONFLICT = None
	if dOriPaths:
		outConflicts = os.path.join(outDir, "%s.agouti_vs_original.compare.txt" %(prefix))
		fCONFLICT = open(outConflicts, 'w')
	with open(outScaffPath, 'w') as fSCAFPATH:
		for i, path in enumerate(agPaths):
			scafName = prefix + "_scaf_%d" %(i)
			fSCAFPATH.write(">%s\n%s\n" %(scafName, ",".join(path)))
			if dOriPaths:
				conflictType = agSUM.check_consistency(dOriPaths, path)
				if conflictType is not None:
					fCONFLICT.write("%s\t%s\t%s\n" %(scafName, conflictType,
													 ",".join(path)))

def summarize_gene_path(dMergedGene2Genes, dMergedGene2Ctgs,
						outDir, prefix):
	'''
		output what each merged gene
		was made of
	'''
	outGenePath = os.path.join(outDir, "%s.agouti.gene_path.txt" %(prefix))
	with open(outGenePath, 'w') as fGENEPATH:
		for k, v in sorted(dMergedGene2Ctgs.iteritems()):
			fGENEPATH.write(">%s\nCONTIGPATH\t%s\nGENEPATH\t%s\n"
							%(k, ','.join(v), ','.join(dMergedGene2Genes[k])))

def output_gff(dGeneModels, dMergedGene2Ctgs, dMergedGene2Genes,
			   outDir, prefix):
	'''
		output all gene models as GFF
	'''
	outGFF = os.path.join(outDir, "%s.agouti.gff" %(prefix))
	numGenes = 0
	with open(outGFF, 'w') as fOUTGFF:
		fOUTGFF.write("##gff-version3\n")
		fOUTGFF.write("# This output was generated with AGOUTI (version 0.1)\n")
		for k, v in dGeneModels.iteritems():
			for i in range(len(v)):
				geneModel = v[i]
				if geneModel.fake == 1:
					continue
				else:
					numGenes += 1
					fOUTGFF.write("# start gene %s\n" %(geneModel.geneID))
					if geneModel.geneID.startswith("AGOUTI_Merged"):
						fOUTGFF.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;MERGE_FROM_GENES=%s;MERGE_FROM_CONTIGS=%s\n" %(k, geneModel.program, geneModel.geneStart, geneModel.geneStop, geneModel.strand, geneModel.geneID, ",".join(dMergedGene2Genes[geneModel.geneID]), ",".join(dMergedGene2Ctgs[geneModel.geneID])))
					else:
						fOUTGFF.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s\n"
									  %(k, geneModel.program, geneModel.geneStart,
										geneModel.geneStop, geneModel.strand,
										geneModel.geneID))
					for j in range(0, len(geneModel.lcds), 2):
						cdsStart = geneModel.lcds[j]
						cdsStop = geneModel.lcds[j+1]
						fOUTGFF.write("%s\t%s\tExon\t%d\t%d\t.\t%s\t.\tID=%s.exon\n"
									  %(k, geneModel.program, cdsStart, cdsStop,
										geneModel.strand, geneModel.geneID))
					fOUTGFF.write("# end gene %s\n" %(geneModel.geneID))
					fOUTGFF.write("###\n")
	return numGenes

def get_gene_index(geneModels, curGeneID, debug=0, reverse=False):
	currentGeneIndex = -1
	if reverse:
		currentGeneIndex = [x.geneID for x in geneModels[::-1]].index(curGeneID)
		if debug:
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\treversed gene models - %s"
										 %(" ".join([x.geneID for x in geneModels[::-1]])))
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\tcurGene - %s - index - %d"
										 %(curGeneID, currentGeneIndex))
	else:
		print ">>>> get_gene_index", " ".join([x.geneID for x in geneModels]), "curGeneID", curGeneID
		currentGeneIndex = [x.geneID for x in geneModels].index(curGeneID)
		if debug:
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\tgene models - %s"
										 %(" ".join([x.geneID for x in geneModels[::-1]])))
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\tcurGene - %s - index - %d"
										 %(curGeneID, currentGeneIndex))
	toLeft = currentGeneIndex
	toRight = len(geneModels) - 1 - currentGeneIndex
	return toLeft, toRight

def merge_gene_model(currentGene, nextGene, scafName,
					 numMergedGene, currentOffset, nextOffset,
					 debug=0):
	mergedGene = agGFF.AGOUTI_GFF()
	mergedGene.geneID = "AGOUTI_Merged_gene_%d" %(numMergedGene)
	mergedGene.geneStart = currentGene.geneStart + currentOffset
	mergedGene.geneStop = nextGene.geneStop + nextOffset
	mergedGene.program = "AGOUTI"
	mergedGene.ctgID = scafName
	if currentGene.fake:
		if nextGene.fake:
			mergedGene.strand = '.'
		else:
			mergedGene.strand = nextGene.strand
	else:
		if nextGene.fake:
			mergedGene.strand = currentGene.strand
		else:
			if currentGene.strand == nextGene.strand:
				mergedGene.strand = currentGene.strand
			else:
				if currentGene.merge:
					mergedGene.strand = nextGene.strand
				else:
					mergedGene.strand = '.'
#	if currentGene.strand == nextGene.strand:
#		mergedGene.strand = currentGene.strand
#	else:
#		mergedGene.strand = '.'
	if debug:
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tmerge cur GeneID - %s - %s"
									 %(currentGene.geneID, str(currentGene.lcds)))
	for i in range(len(currentGene.lcds)):
		currentGene.lcds[i] += currentOffset
	for i in range(len(nextGene.lcds)):
		nextGene.lcds[i] += nextOffset
	mergedGene.lcds = currentGene.lcds + nextGene.lcds
	mergedGene.merge = 1
	if debug:
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tmerge next GeneID - %s - %s"
									 %(nextGene.geneID, str(nextGene.lcds)))
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tmerged GeneID - %s - %s"
									 %(mergedGene.geneID, str(mergedGene.lcds)))
	#print ">>>> merge", "next LCDS after", nextGene.lcds
	#print ">>>> merge", "finish merging", mergedGene.lcds
	return mergedGene

def update_gene_model(geneModels, updatedGeneModels,
					  newScafID, offset, excludeIDs,
					  debug=0, mergedGene=None):

	if debug:
		agUPDATEDebug.debugger.debug("UPDATE_GENE_MODEL\t\texcludes to update - %s"
									 %(str(excludeIDs)))
	tmpGeneModel = agGFF.AGOUTI_GFF()
	for i in range(len(geneModels)):
		tmpGeneModel = geneModels[i]
		if tmpGeneModel.geneID not in excludeIDs:
			tmpGeneModel = geneModels[i]
			tmpGeneModel.ctgID = newScafID
			tmpGeneModel.geneStart += offset
			tmpGeneModel.geneStop += offset
			for j in range(len(tmpGeneModel.lcds)):
				tmpGeneModel.lcds[j] += offset
			updatedGeneModels.append(tmpGeneModel)
			#print ">>>>", "update", tmpGeneModel.geneID, tmpGeneModel.lcds
	if mergedGene is not None:
		#print ">>>> mergedGene to update", mergedGene.geneID
		updatedGeneModels.append(mergedGene)
		if debug:
			agUPDATEDebug.debugger.debug("UPDATE_GENE_MODEL\t\tUpdate gene - %s - %s"
										 %(mergedGene.geneID, str(mergedGene.lcds)))
		#print ">>>>", "update", mergedGene.geneID, mergedGene.lcds
	updatedGeneIDs = []
	for i in range(len(updatedGeneModels)):
		updatedGeneIDs.append(updatedGeneModels[i].geneID)
	if debug:
		agUPDATEDebug.debugger.debug("UPDATE_GENE_MODEL\t\tupdated GeneIDs - %s"
									 %(updatedGeneIDs))
	return updatedGeneModels, updatedGeneIDs

def reverse_gene_models(geneModels, lenCtg, debug=0, excludeIDs=[]):
	'''
		reverse gene models
		call reverse_gene_model() function
	'''
	if debug:
		agUPDATEDebug.debugger.debug("REV_GENE_MODELS\t\texcludes to reverse - %s"
									 %(str(excludeIDs)))
	for i in xrange(len(geneModels)):
		if geneModels[i].geneID not in excludeIDs:
			geneModels[i] = reverse_gene_model(geneModels[i], lenCtg, debug)
	return geneModels

def exclude_gene_model(geneModels, exclude):
	tmp = []
	for i in xrange(len(geneModels)):
		if geneModels[i].geneID != exclude:
			tmp.append(geneModels[i])
	return tmp

def reverse_gene_model(geneModel, lenCtg, debug=0):
	'''
		reverse each gene, all CDS within the gene
	'''
	tmpGeneStop =  lenCtg-geneModel.geneStart+1
	tmpGeneStart = lenCtg-geneModel.geneStop+1
	geneModel.geneStart = tmpGeneStart
	geneModel.geneStop = tmpGeneStop
	lcds = geneModel.lcds
	reverseLCDs = []
	if debug:
		agUPDATEDebug.debugger.debug("REV_GENE_MODEL\t\tbefore reverse - %s - %s"
										%(geneModel.geneID, str(lcds)))
	for j in range(len(lcds)-1, -1, -2):
		reverseLCDs += [lenCtg-lcds[j]+1, lenCtg-lcds[j-1]+1]
	if debug:
		agUPDATEDebug.debugger.debug("REV_GENE_MODEL\t\tafter reverse - %s - %s"
										%(geneModel.geneID, str(reverseLCDs)))
	geneModel.lcds = reverseLCDs
	if geneModel.strand == '+':
		geneModel.strand = '-'
	else:
		geneModel.strand = '+'
	return geneModel

def ctgpair2genepair(dCtgPair2GenePair, curCtg, nextCtg):
	currentGene = None
	nextGene = None
	if (curCtg, nextCtg) in dCtgPair2GenePair:
		genePair = dCtgPair2GenePair[curCtg, nextCtg]
		currentGene = genePair[0]
		nextGene = genePair[1]
	elif (nextCtg, curCtg) in dCtgPair2GenePair:
		genePair = dCtgPair2GenePair[nextCtg, curCtg]
		currentGene = genePair[1]
		nextGene = genePair[0]
	return currentGene, nextGene

def get_orientation_counts(curVertex, nextVertex, edgeSenseDict, debug=0):
	''' get the counts of each orientation '''
	senseList = []
	if (curVertex, nextVertex) in edgeSenseDict:
		senseList = edgeSenseDict[curVertex, nextVertex]
		FR = senseList.count(("+", "-"))
		RF = senseList.count(("-", "+"))
	elif (nextVertex, curVertex) in edgeSenseDict:
		senseList = edgeSenseDict[nextVertex, curVertex]
		# flip because the from and end vertices fliped
		FR = senseList.count(("-", "+"))
		RF = senseList.count(("+", "-"))
#	else:
#		if debug:
#			agUPDATEProgress.debugger.debug("missing orientation for these two nodes: %d, %d"
#											%(curVertex, nextVertex)
#		return -1, -1, -1, -1
	FF = senseList.count(("+", "+"))
	RR = senseList.count(("-", "-"))
	return FR, FF, RR, RF

def decide_orientation(FR, FF, RR, RF):
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
