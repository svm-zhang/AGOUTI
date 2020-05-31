import os
import sys
import collections

from lib import agouti_log as agLOG
from lib import agouti_gff as agGFF
from src import agouti_sequence as agSeq
from src import agouti_path as agPATH

def agouti_update(agoutiPaths, dSeqs, vertex2Name,
				  dSenses, dGFFs,
				  dCtgPair2GenePair, outDir, prefix,
				  nFills=1000, debug=0, no_update_gff=0):

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

	if not no_update_gff:
		agUPDATEProgress.logger.info("[BEGIN] Updating gene models")

	outFasta = os.path.join(outDir, "%s.agouti.fasta" %(prefix))
	fFASTA = open(outFasta, 'w')
	dUpdateGFFs = collections.defaultdict(list)
	dMergedGene2Ctgs = collections.defaultdict(list)
	dMergedGene2Genes = collections.defaultdict(list)
	scafPaths = []
	numMergedGene = 0
	nCtgScaffolded = 0
	scaffoldedCtgs = {}
	seqLens = []
	dScafGaps = {}
	dScafStats = {}
	scafID = 0
	mergedGenes = []
	for i in range(len(agoutiPaths)):
		path = agoutiPaths[i]
		scafID += 1
		scafName = prefix + "_scaf_%d" %(scafID)
		dScafStats[scafName] = 0
		dScafGaps[scafName] = []

		curVertex = path[0]
		sequence = dSeqs[curVertex]
		curSense = "+"
		curCtg = vertex2Name[curVertex]
		preCtg = ""
		scafPath = [curVertex]
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
			nextCtg = vertex2Name[nextVertex]

			if preCtg == "":
				if debug:
					agUPDATEDebug.debugger.debug("UPDATE_MAIN\t>scaf_%d - path - %s"
												 %(scafID,
												  str([vertex2Name[vertex] for vertex in path])))
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tcurVertex - %d - %s - nextVertex - %d - %s"
											 %(curVertex, curCtg, nextVertex, nextCtg))

			if not no_update_gff:
				#curGene, nextGene = ctgpair2genepair(dCtgPair2GenePair, curCtg, nextCtg)
				curGene, nextGene = ctgpair2genepair(dCtgPair2GenePair, curVertex, nextVertex)
				#!!! I should not break here, should continue#
				if curGene is None and nextGene is None:
					agUPDATEProgress.logger.error("%s - %s found no gene models joining them"
										   %(curCtg, nextCtg))
					agUPDATEProgress.logger.error("This is NOT EXPECTED, REPORT!")
					sys.exit(1)
				curGeneID = curGene.geneID
				excludeGeneIDs = [preGeneID] + [curGeneID]
				if debug:
					agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tpreGene - %s - curGene - %s - nextGene - %s"
												 %(preGeneID, curGene.geneID, nextGene.geneID))

			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafName - %s" %(scafName))
			FR, FF, RR, RF = get_orientation_counts(curVertex, nextVertex, dSenses)
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tcurSense=%s" %(curSense))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tFR=%d - FF=%d - RF=%d - RR=%d"
											 %(FR, FF, RF, RR))
			if curSense == "-":
				temp1 = FR
				temp2 = FF
				FR = RR
				FF = RF
				RR = temp1
				RF = temp2
			orientation = decide_orientation(FR, FF, RR, RF)

			gapStart = gapStop + len(dSeqs[curVertex])
			gapStop = gapStart + nFills - 1
			dScafGaps[scafName].append((gapStart+1, gapStop+1))
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tcurSense=%s" %(curSense))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tFR=%d - FF=%d - RF=%d - RR=%d"
											 %(FR, FF, RF, RR))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\toffset - %d - curCtgLen - %d"
											 %(offset, len(dSeqs[curVertex])))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tgapstart - %d - gapstop - %d"
											 %(gapStart, gapStop+1))
			valid = 0
			if orientation == "FR":
				if not no_update_gff:
					if curGeneID !=  preGeneID:
						numMergedGene += 1
						mergedGene = merge_gene_model(curGene, nextGene, scafName,
													  numMergedGene, offset, gapStop,
													  debug)
						dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
						#if curGene.geneStop != 0:
						#	dMergedGene2Genes[mergedGene.geneID] += [curGeneID]
						#if nextGene.geneStop != 0:
						#	dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
						if mergedGene.geneStop != 0:
							dMergedGene2Genes[mergedGene.geneID] += [curGeneID, nextGene.geneID]
						dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg], dUpdateGFFs[scafName],
																				  scafName, offset, excludeGeneIDs,
																				  debug, mergedGene)
					else:
						mergedGene = merge_gene_model(preMergedGene, nextGene, scafName,
													  numMergedGene, 0, gapStop, debug)
						dMergedGene2Ctgs[mergedGene.geneID] += [nextCtg]
						if nextGene.geneStop != 0:
							dMergedGene2Genes[mergedGene.geneID] += [nextGene.geneID]
						indexMerged = updatedGeneIDs.index(mergedGene.geneID)
						dUpdateGFFs[scafName][indexMerged] = mergedGene
					preMergedGene = mergedGene
				sequence += 'N'*nFills + dSeqs[nextVertex]
				scafPath += [nextVertex]
				curSense = "+"
			elif orientation == "FF":
				if not no_update_gff:
					#nextGene = reverse_gene_model(nextGene, len(dSeqs[nextVertex]), debug)
					dGFFs[nextCtg] = reverse_gene_models(dGFFs[nextCtg], len(dSeqs[nextVertex]), debug)
					if curGeneID !=  preGeneID:
						numMergedGene += 1
						mergedGene = merge_gene_model(curGene, nextGene, scafName,
													  numMergedGene, offset, gapStop,
													  debug)
						dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
						if mergedGene.geneStop != 0:
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
				sequence += 'N'*nFills + agSeq.rc_seq(dSeqs[nextVertex])
				scafPath += [-1*nextVertex]
				curSense = "-"
			elif orientation == "RR":
				if not no_update_gff:
					if curGene.geneID != preGeneID:
						dGFFs[curCtg] = reverse_gene_models(dGFFs[curCtg], len(dSeqs[curVertex]), debug)
						numMergedGene += 1
						mergedGene = merge_gene_model(curGene, nextGene, scafName,
													  numMergedGene, offset, gapStop,
													  debug)
						dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
						if mergedGene.geneStop != 0:
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
					preMergedGene = mergedGene
				sequence = agSeq.rc_seq(sequence) + \
						   'N'*nFills + dSeqs[nextVertex]
				scafPath[-1] = -1*scafPath[-1]
				scafPath += [nextVertex]
				curSense = "+"
			elif orientation == "RF":
				if not no_update_gff:
					dGFFs[nextCtg] = reverse_gene_models(dGFFs[nextCtg], len(dSeqs[nextVertex]), debug)
					if curGene.geneID != preGeneID:
						dGFFs[curCtg] = reverse_gene_models(dGFFs[curCtg], len(dSeqs[curVertex]), debug)
						#nextGene = reverse_gene_model(nextGene, len(dSeqs[nextVertex]), debug)
						numMergedGene += 1
						mergedGene = merge_gene_model(curGene, nextGene, scafName,
													  numMergedGene, offset, gapStop,
													  debug)
						dMergedGene2Ctgs[mergedGene.geneID] += [curCtg, nextCtg]
						if mergedGene.geneStop != 0:
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
					preMergedGene = mergedGene
				sequence = agSeq.rc_seq(sequence) + \
						   'N'*nFills + \
						   agSeq.rc_seq(dSeqs[nextVertex])
				scafPath[-1] = -1*scafPath[-1]
				scafPath += [-1*nextVertex]
				curSense = "-"
			if debug:
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafPath in vertices updates- %s"
											 %(str(scafPath)))
				agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tdMergedGene2Gene - %s"
											 %(str(dMergedGene2Genes[mergedGene.geneID])))
			if not no_update_gff:
				mergedGenesPerPath.append(mergedGene.geneID)
				preGeneID = nextGene.geneID
			offset = gapStop
			preCtg = curCtg
			curVertex = nextVertex
			curCtg = vertex2Name[curVertex]

		for i in range(len(scafPath)):
			v = scafPath[i]
			if v < 0:
				scafPath[i] = "-"+vertex2Name[-1*v]
			else:
				scafPath[i] = vertex2Name[v]
		scafPaths += [scafPath]
		if debug:
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafPath in human-readable updates- %s"
										 %(str(scafPath)))
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tappend last curCtg - %s" %(curCtg))
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tscafPath - %s"
										 %(str(scafPath)))
		if not no_update_gff:
			excludeGeneIDs = [preGeneID]
			mergedGenes.append(mergedGenesPerPath)
			dUpdateGFFs[scafName], updatedGeneIDs = update_gene_model(dGFFs[curCtg],
																	  dUpdateGFFs[scafName],
																	  scafName, offset,
																	  excludeGeneIDs, debug)
		fFASTA.write(">%s |%dbp |%s\n%s\n"
						%(scafName, len(sequence), ",".join(scafPath), sequence))
		dScafStats[scafName] = len(sequence)
		seqLens.append(len(sequence))
		#agPaths.append(scafPath)
		nCtgScaffolded += len(scafPath)
		scaffoldedCtgs.update(dict((contig, 1) for contig in scafPath))
		if debug:
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t\tmergedGenesPerPath - %s"
										 %(str(mergedGenesPerPath)))
			agUPDATEDebug.debugger.debug("UPDATE_MAIN\t-------------------------------------")

	agPATH.report_scaffold_path(scafPaths, vertex2Name, outDir, prefix)

	# other contigs need to be output
	agUPDATEProgress.logger.info("Finalizing sequences")
	for vertex in dSeqs:
		if vertex2Name[vertex] in scaffoldedCtgs or "-"+vertex2Name[vertex] in scaffoldedCtgs:
			continue
		fFASTA.write(">%s\n%s\n" % (vertex2Name[vertex], dSeqs[vertex]))
		dScafStats[vertex2Name[vertex]] = len(dSeqs[vertex])
		seqLens.append(len(dSeqs[vertex]))
	fFASTA.close()
	n50 = agSeq.get_assembly_NXX(seqLens)

	agUPDATEProgress.logger.info("Outputting updated Gene Moddels")
	for vertex in dSeqs:
		if vertex2Name[vertex] in scaffoldedCtgs:
			if vertex2Name[vertex] in dGFFs:
				del dGFFs[vertex2Name[vertex]]
	if not no_update_gff:
		dFinalGFFs = dict(dGFFs, **dUpdateGFFs)
		numGenes = output_gff(dFinalGFFs, dMergedGene2Ctgs, dMergedGene2Genes,
							  dScafStats, dScafGaps, outDir, prefix)
		agUPDATEProgress.logger.info("Summarizing AGOUTI gene paths")
		summarize_gene_path(dMergedGene2Genes, dMergedGene2Ctgs,
							outDir, prefix)

	agUPDATEProgress.logger.info("-----------Summary-----------")
	agUPDATEProgress.logger.info("number of contigs scaffoled: %d" %(nCtgScaffolded))
	agUPDATEProgress.logger.info("number of scaffolds: %d" %(scafID))
	agUPDATEProgress.logger.info("number of contigs in the final assembly: %d" %(len(seqLens)))
	agUPDATEProgress.logger.info("Final assembly N50: %d" %(n50))
	if not no_update_gff:
		agUPDATEProgress.logger.info("Final number of genes: %d" %(numGenes))
	agUPDATEProgress.logger.info("Succeeded")

def recover_original_scaffold():
	'''
		not finish
	'''
	if dOriPaths is None:
		for vertex in sorted(dSeqs):
			if vertex2Name[vertex] not in scaffoldedCtgs:
				fFASTA.write(">%s\n%s\n" % (vertex2Name[vertex], dSeqs[vertex]))
				seqLens.append(len(dSeqs[vertex]))
	else:
		for _, oriPath in dOriPaths.items():
			if len(oriPath) == 1:
				index = vertex2Name.index(oriPath[0])
				fFASTA.write(">%s\n%s\n" % (vertex2Name[index], dSeqs[index]))
				seqLens.append(len(dSeqs[index]))
				continue
			print(oriPath)
			untouchedCtgs = [k for k in oriPath if k not in scaffoldedCtgs]
			print(untouchedCtgs)
			sequence = ""
			preCtg = untouchedCtgs[0]
			preIndex = vertex2Name.index(preCtg)
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
				preIndex = vertex2Name.index(preCtg)
				curIndex = vertex2Name.index(curCtg)
				sequence += gapLen * 'N' + \
							dSeqs[curIndex]
				print("preCtg", preCtg, "curCtg", curCtg)
				print(len(sequence))
				preCtg = curCtg
			if sequence:
				fFASTA.write(">%s\n%s\n" % ("haha_need_to_fix", sequence))
				seqLens.append(len(sequence))
			sys.exit()

def summarize_gene_path(dMergedGene2Genes, dMergedGene2Ctgs,
						outDir, prefix):
	'''
		output what each merged gene
		was made of
	'''
	outGenePath = os.path.join(outDir, "%s.agouti.gene_path.txt" %(prefix))
	with open(outGenePath, 'w') as fGENEPATH:
		for k, v in sorted(dMergedGene2Ctgs.items()):
			# only actually merged gene output
			if dMergedGene2Genes[k]:
				fGENEPATH.write(">%s\nCONTIGPATH\t%s\nGENEPATH\t%s\n"
								%(k, ','.join(v), ','.join(dMergedGene2Genes[k])))

def output_gff(dGeneModels, dMergedGene2Ctgs, dMergedGene2Genes,
			   dScafStats, dScafGaps, outDir, prefix):
	'''
		output all gene models as GFF
	'''
	outGFF = os.path.join(outDir, "%s.agouti.gff" %(prefix))
	numGenes = 0
	dSeen = {}
	with open(outGFF, 'w') as fOUTGFF:
		fOUTGFF.write("##gff-version3\n")
		fOUTGFF.write("# This output was generated with AGOUTI (version 0.3.1)\n")
		for k, v in dGeneModels.items():
			if k not in dScafStats:
				continue
			if k not in dSeen:
				fOUTGFF.write("%s\tAGOUTI\tscaffold\t1\t%d\t.\t.\t.\tID=%s\n"
							  %(k, dScafStats[k], k))
				if k in dScafGaps:
					gaps = dScafGaps[k]
					for i in range(len(gaps)):
						fOUTGFF.write("%s\tAGOUTI\tgap\t%d\t%d\t.\t.\t.\tParent=%s\n"
									  %(k, gaps[i][0], gaps[i][1], k))
			tmpV = sorted([(i.geneStart, i.geneStop, i) for i in v])
			for i in range(len(tmpV)):
				geneModel = tmpV[i][2]
				# no output of fake gene created
				if geneModel.fake == 1:
					continue
				# no output of gene with zero intervals
				# this is designed for recovering original
				# assembly without any gene models
				elif geneModel.geneStart == 0 and geneModel.geneStop == 0:
					continue
				else:
					numGenes += 1
					fOUTGFF.write("# start gene %s\n" %(geneModel.geneID))
					if geneModel.geneID.startswith("agMerge"):
						fOUTGFF.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;MERGE_FROM_GENES=%s;MERGE_FROM_CONTIGS=%s\n"
									  %(k, geneModel.program, geneModel.geneStart,
										geneModel.geneStop,
										geneModel.strand, geneModel.geneID,
										",".join(dMergedGene2Genes[geneModel.geneID]),
										",".join(dMergedGene2Ctgs[geneModel.geneID])))
					else:
						fOUTGFF.write("%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s\n"
									  %(k, geneModel.program, geneModel.geneStart,
										geneModel.geneStop, geneModel.strand,
										geneModel.geneID))
					if len(geneModel.lcds) >= 2:
						mrnaStart = min(geneModel.lcds[0], geneModel.lcds[-1])
						mrnaStop = max(geneModel.lcds[0], geneModel.lcds[-1])
						mrnaID = geneModel.geneID+".mRNA"
						fOUTGFF.write("%s\t%s\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n"
									  %(k, geneModel.program, mrnaStart,
										mrnaStop, geneModel.strand,
										mrnaID, geneModel.geneID))
					for j in range(0, len(geneModel.lcds), 2):
						cdsStart = geneModel.lcds[j]
						cdsStop = geneModel.lcds[j+1]
						fOUTGFF.write("%s\t%s\texon\t%d\t%d\t.\t%s\t.\tID=%s.exon.%d;Parent=%s\n"
									  %(k, geneModel.program, cdsStart, cdsStop,
									  geneModel.strand, geneModel.geneID, (j/2)+1, mrnaID))
					fOUTGFF.write("# end gene %s\n" %(geneModel.geneID))
					fOUTGFF.write("###\n")
	return numGenes

def get_gene_index(geneModels, curGeneID, debug=0, reverse=False):
	curGeneIndex = -1
	if reverse:
		curGeneIndex = [x.geneID for x in geneModels[::-1]].index(curGeneID)
		if debug:
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\treversed gene models - %s"
										 %(" ".join([x.geneID for x in geneModels[::-1]])))
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\tcurGene - %s - index - %d"
										 %(curGeneID, curGeneIndex))
	else:
		print(">>>> get_gene_index", " ".join([x.geneID for x in geneModels]), "curGeneID", curGeneID)
		curGeneIndex = [x.geneID for x in geneModels].index(curGeneID)
		if debug:
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\tgene models - %s"
										 %(" ".join([x.geneID for x in geneModels[::-1]])))
			agUPDATEDebug.debugger.debug("GET_GENE_INDEX\t\tcurGene - %s - index - %d"
										 %(curGeneID, curGeneIndex))
	toLeft = curGeneIndex
	toRight = len(geneModels) - 1 - curGeneIndex
	return(toLeft, toRight)

def merge_gene_model(curGene, nextGene, scafName,
					 numMergedGene, currentOffset, nextOffset,
					 debug=0):
	'''
		merge two gene models
	'''
	mergedGene = agGFF.AGOUTI_GFF()
	mergedGene.geneID = "agMerge_%d" %(numMergedGene)
	if curGene.geneStop == 0 and nextGene.geneStop == 0:
		mergedGene.geneStart = 0
		mergedGene.geneStop = 0
	elif curGene.geneStop == 0 and nextGene.geneStop > 0:
		mergedGene.geneStart = nextGene.geneStart + nextOffset
		mergedGene.geneStop = nextGene.geneStop + nextOffset
	elif curGene.geneStop != 0 and nextGene.geneStop == 0:
		mergedGene.geneStart = curGene.geneStart + currentOffset
		mergedGene.geneStop = curGene.geneStop + currentOffset
	elif curGene.geneStop != 0 and nextGene.geneStop != 0:
		mergedGene.geneStart = curGene.geneStart + currentOffset
		mergedGene.geneStop = nextGene.geneStop + nextOffset
	#mergedGene.geneStart = curGene.geneStart + currentOffset
	#mergedGene.geneStop = nextGene.geneStop + nextOffset
	mergedGene.program = "AGOUTI"
	mergedGene.ctgID = scafName
	if curGene.fake:
		if nextGene.fake:
			mergedGene.strand = '.'
		else:
			mergedGene.strand = nextGene.strand
	else:
		if nextGene.fake:
			mergedGene.strand = curGene.strand
		else:
			if curGene.strand == nextGene.strand:
				mergedGene.strand = curGene.strand
			else:
				if curGene.merge:
					mergedGene.strand = nextGene.strand
				else:
					mergedGene.strand = '.'
	#if debug:
	for i in range(len(curGene.lcds)):
		curGene.lcds[i] += currentOffset
	for i in range(len(nextGene.lcds)):
		nextGene.lcds[i] += nextOffset
	mergedGene.lcds = curGene.lcds + nextGene.lcds
	mergedGene.merge = 1
	if debug:
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tcuroffset=%d - nextoffset=%d"
									 %(currentOffset, nextOffset))
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tcur GeneID - %s - %d - %d - %s"
									 %(curGene.geneID, curGene.geneStart+currentOffset,
									 curGene.geneStop+currentOffset, str(curGene.lcds)))
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tnext GeneID - %s - %d - %d - %s"
									 %(nextGene.geneID, nextGene.geneStart+nextOffset,
									 nextGene.geneStop+nextOffset, str(nextGene.lcds)))
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tmerged GeneID - %s - %d - %d - %s"
									 %(mergedGene.geneID, mergedGene.geneStart,
									 mergedGene.geneStop, str(mergedGene.lcds)))
		agUPDATEDebug.debugger.debug("MERGE_GENE_MODEL\t\tmergeGeneStart - %d mergedGeneStop - %d"
									 %(mergedGene.geneStart, mergedGene.geneStop))
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
	for i in range(len(geneModels)):
		if geneModels[i].geneID not in excludeIDs:
			geneModels[i] = reverse_gene_model(geneModels[i], lenCtg, debug)
	return geneModels

def reverse_gene_model(geneModel, lenCtg, debug=0):
	'''
		reverse each gene, all CDS within the gene
	'''
	if debug:
		agUPDATEDebug.debugger.debug("REV_GENE_MODEL\t\tbefore reverse - %s - %d - %d - %s"
										%(geneModel.geneID, geneModel.geneStart,
										  geneModel.geneStop, str(geneModel.lcds)))
	tmpGeneStop =  lenCtg-geneModel.geneStart+1
	tmpGeneStart = lenCtg-geneModel.geneStop+1
	geneModel.geneStart = tmpGeneStart
	geneModel.geneStop = tmpGeneStop
	lcds = geneModel.lcds
	reverseLCDs = []
	for j in range(len(lcds)-1, -1, -2):
		reverseLCDs += [lenCtg-lcds[j]+1, lenCtg-lcds[j-1]+1]
	if debug:
		agUPDATEDebug.debugger.debug("REV_GENE_MODEL\t\tafter reverse - %s - %d - %d - %s"
										%(geneModel.geneID, geneModel.geneStart,
										  geneModel.geneStop, str(reverseLCDs)))
	geneModel.lcds = reverseLCDs
	if geneModel.strand == '+':
		geneModel.strand = '-'
	else:
		geneModel.strand = '+'
	return geneModel

def exclude_gene_model(geneModels, exclude):
	tmp = []
	for i in range(len(geneModels)):
		if geneModels[i].geneID != exclude:
			tmp.append(geneModels[i])
	return tmp

def ctgpair2genepair(dCtgPair2GenePair, curCtg, nextCtg):
	curGene = None
	nextGene = None
	if (curCtg, nextCtg) in dCtgPair2GenePair:
		genePair = dCtgPair2GenePair[curCtg, nextCtg]
		curGene = genePair[0]
		nextGene = genePair[1]
	elif (nextCtg, curCtg) in dCtgPair2GenePair:
		genePair = dCtgPair2GenePair[nextCtg, curCtg]
		curGene = genePair[1]
		nextGene = genePair[0]
	return curGene, nextGene

def get_orientation_counts(curVertex, nextVertex, dSenses, debug=0):
	''' get the counts of each orientation '''
	senseList = []
	if (curVertex, nextVertex) in dSenses:
		senseList = dSenses[curVertex, nextVertex]
		FR = senseList.count(("+", "-"))
		RF = senseList.count(("-", "+"))
	elif (nextVertex, curVertex) in dSenses:
		senseList = dSenses[nextVertex, curVertex]
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
