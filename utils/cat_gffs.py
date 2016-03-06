import os
import sys

def cat(gffFiles, outFile):
	fOUT = open(outFile, 'w')
	nGenes = 1
	for (i, gff) in enumerate(gffFiles):
		print os.path.realpath(gff)
		with open(gff, 'r') as fGFF:
			if i != 0:
				for i in range(12):
					fGFF.readline()
			n = 1
			for line in fGFF:
				tmpLine = line.strip().split("\t")
				if line.startswith("# start"):
					fOUT.write("# start gene g%d\n" %(nGenes))
				elif line.startswith("# end"):
					fOUT.write("# end gene g%d\n" %(nGenes))
					nGenes += 1
				elif len(tmpLine) == 9:
					if tmpLine[2] == "gene":
						fOUT.write("%s\tID=g%d\n" %("\t".join(tmpLine[:-1]), nGenes))
					elif tmpLine[2] == "transcript":
						fOUT.write("%s\tID=g%d.t1;Parent=g%d\n" %("\t".join(tmpLine[:-1]), nGenes, nGenes))
					elif tmpLine[2] == "transcription_start_site":
						fOUT.write("%s\tParent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes))
					elif tmpLine[2] == "start_codon":
						fOUT.write("%s\tParent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes))
					elif tmpLine[2] == "intron":
						fOUT.write("%s\tParent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes))
					elif tmpLine[2] == "CDS":
						fOUT.write("%s\tID=g%d.t1.cds;Parent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes, nGenes))
					elif tmpLine[2] == "exon":
						fOUT.write("%s\tParent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes))
					elif tmpLine[2] == "stop_codon":
						fOUT.write("%s\tParent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes))
					elif tmpLine[2] == "transcription_end_site":
						fOUT.write("%s\tParent=g%d.t1\n" %("\t".join(tmpLine[:-1]), nGenes))
				else:
					fOUT.write(line)
	print nGenes - 1

def main():
	if len(sys.argv) < 3:
		print "python cat_gffs.py *.gff cat.gff"
		sys.exit(1)
	gffFiles = sys.argv[1:-1]
	outFile = sys.argv[-1]
	cat(sys.argv[1:-1], sys.argv[-1])

if __name__ == "__main__":
	main()
