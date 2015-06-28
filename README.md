Annotated Genome Optimization Using Transcriptome Information (AGOUTI)
==========================
AGOUTI v0.1

Description
=======================
AGOUTI uses RNAseq reads to guide genome scaffolding and improve gene annotation.

```
AGOUTI works in the following steps:
	i) read RNAseq mapping results and get reads pairs uniquely mapped to different contigs
	ii) filter reads pairs based on their mapped positions and orientations
	iii) scaffolding based on reads pairs passing step ii filtration, this step uses a modified version of RNAPATH
	iv) updating assembly and gene annotations
```

Requirement
==========================
1. python 2.7 or above
2. SAMtools

Running AGOUTI
==========================
```
python agouti.py -h
usage: agouti.py [-h] -contig FILE -bam FILE -gff FILE -out DIR [-mnl INT] [-p STR]

Welcome to AGOUTI!

optional arguments:
	-h, --help    show this help message and exit
	-contig FILE  specify the initial assembly in FASTA format
	-bam FILE     specify the RNA-seq mapping results in BAM format
	-gff FILE     specify the predicted gene model in GFF format
	-out DIR      specify the directory to store output files
	-mnl INT      minimum number of reads supporting a link between a contig pair
	-p STR        specify the output prefix
```

Example
==========================
```
samtools view test.bam | \
python agouti.py | \
-contig test.fasta | \
-bam - | \
-gff test.gff | \
-out ./test | \
-p agouti > stdout
```

Output
==========================
```
Under test folder, you will find:
	test.agouti.fasta: scaffolded assembly in FASTA format
	test.agouti.gff3: updated gene annotations in GFF3 format
	test.scaff.paths: details on scaffolding paths
	test.join_pairs: reads pairs used for scaffolding
```
