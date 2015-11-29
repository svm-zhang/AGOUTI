import sys
import math

def read_original_path(oriScafPathFile):
	dOriPaths = {}
	dOriGaps = {}
	with open(oriScafPathFile, 'r') as fIN:
		oriPath = []
		contigID = ""
		for line in fIN:
			if line.startswith('>'):
				if oriPath:
					dOriPaths[oriScafID] = oriPath
				oriScafID = line.strip()[1:]
				oriPath = []
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
				if oriPath:
					oriPath.append(nexCtg)
				else:
					oriPath += [curCtg, nexCtg]
		if oriPath:
			dOriPaths[oriScafID] = oriPath

#	print "Number of scaffolds in the original assembly: %d" %(len(dOriPaths))
	return dOriPaths, dOriGaps

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

