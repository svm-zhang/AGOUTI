# **AGOUTI**: Annotated Genome Optimization Using Transcriptome Information

## Overview
AGOUTI(v0.2) uses RNAseq reads to guide genome scaffolding and improve gene annotation.

```
AGOUTI works in the following steps:
	i) read RNAseq mapping results and get reads pairs uniquely mapped to different contigs
	ii) filter reads pairs based on their mapped positions and orientations
	iii) scaffolding based on reads pairs passing step ii filtration, this step uses a modified version of RNAPATH
	iv) updating assembly and gene annotations
```

## Features

1. Scaffold hundreds to thousands of contigs, yielding more contiguous assemblies;
2. Reduce the number of gene models and update them simultaneously;
3. Record any inconsistencies with the original (input) scaffolding results;
4. Support break-and-continue feature such that some time-consuming steps can be skipped if the previous run is successful;
5. Output graph in dot format for visualization

## Obtain AGOUTI

To download AGOUTI, please use git to download the most recent development version through:

    git clone git@github.com:svm-zhang/AGOUTI.git

AGOUTI is designed to be light-weight by only requiring SAMtools and python 2.7 or above.

To get the current version of the download, simply do

    python agouti.py -v

You should see, for example, something like this:

    AGOUTI v0.2

In any case, please use this version as your reference.

## Command-line interface

To get started with AGOUTI, please use `-h` or `--help` for a list of supported options.

    python agouti.py -h

```
usage: agouti.py [-h] -contig FILE -bam FILE -gff FILE -out DIR [-mnl INT] [-p STR]

Welcome to AGOUTI!

optional arguments:
	-h, --help                  show this help message and exit
	-contig             FILE    specify the initial assembly in FASTA format
	-bam                FILE    specify the RNA-seq mapping results in BAM format
	-gff                FILE    specify the predicted gene model in GFF format
    -algorithm          STR     specify the scaffolding algorith: gene or weight priority [gene]
	-outdir             DIR     specify the directory to store output files
	-p                  STR     specify the output prefix [agouti]
	-mnl                INT     minimum number of supporting joining-pairs [5]
    -nN                 INT     number of Ns put in between a pair of contigs [1000]
    -minMapQ            INT     minimum mapping quality to use [5]
    -minFracOvl         FLOAT   minimum fraction of alignment to use [0.0]
    -maxFracMismatch    FLOAT   maximum fraction of mismatch per alignment [1.0]
    -debug                      Output extra info for debug
    -overwrite                  specify to overwirte all results from previous run
    -v, --version               show program's version number and exit
```

## Get Started

In it simplest usage, AGOUTI takes three inputs: an initial genome assembly in FASTA format, paired-end RNA-seq reads mapped against the assembly in BAM format, and gene predictions from the initial assembly in GFF3 format. For instance:

    samtools view example.bam | \
    python agouti.py -contig example.fasta | \
    -bam - | \
    -gff example.gff | \
    -outdir ./example | \
    > ./example.stdout

This will produce a scaffoled assembly in FASTA format, and a updated gene models in GFF3 format. All files (including the intermediate files) will be stored under a directory specified by `-outdir`, "example" in this case.

**Note** AGOUTI does not need header information from the BAM file. So `-h` needs not to be specified.

## Prepare Inputs

Assuming you have a dataset of paired-end RNA-seq reads, `example.1.fq` and `example.2.fq`, and an initial assembly generated from an assembler of your favorite, `example.fasta`. You will first need to map the RNA-seq data against the assembly using a short-reads mapper, such as BWA. For example,

    bwa index example.fasta
    bwa mem example.fasta example.1.fq example.2.fq | samtools view -h - > example.bam

At the end of reads mapping, you will have the mapping results in BAM format.

To run AGOUTI, you will also need a set of gene models predicted from the assembly. For instance,

    Augustus --AUGUSTUS_CONFIG_PATH=[path to augustus config file] -gff3=on --species=[your sepcies] example.fasta > example.gff

At the end of gene prediction, you will now have a set of gene models predicted from the assembly. You can choose any * ab inito * gene predictor as long as it spits out the models in GFF3 format. Specifically, AGOUTI looks for the following information in a GFF3 file.

* a line annotated as `gene`
    * contig ID
    * gene ID
    * start and stop positions of the gene
    * strand
* a line annotated as `CDS`
    * start and stop positions of each coding frame

AGOUTI will not issue any complaints if your gene prediction have these information.

**And that's it! You now have all the inputs required AGOUTI.**

## Understand Outputs

AGOUTI outputs its results to a base directory specified by `-outdir`. Under the base director, there are several sub-folders created, each corresponding to a step built in AGOUTI. As of current version (v0.2), a run of AGOUTI using the command-line setting demonstrated in **Getting Started** will generated a structured output like the following:

* drw-r--r-- example
    * drw-r--r--    agouti_seq
    * drw-r--r--    agouti_gff
    * drw-r--r--    agouti_join_pairs
    * drw-r--r--    scaffolding
    * drw-r--r--    agouti_update
    * -rw-r--r--    agouti.main.log

Each subfolder includes three types of file:

1. general progress meter info
2. debug info
3. intermediate outputs

To get a file with debug info you will need to specify `-debug`. An intermediate file can have all the joining-pairs, the denoised set of joining-pairs, the graph in DOT format, etc. Some intermediate files are important to support the break-and-continue feature, e.g. the file with the denoised-set of joining-pairs (see below for more details).

The `agouti.main.log` is prefixed with the string specified by `-p`, so do all the other files generated by AGOUTI. The sequence ID in the final assembly will also be as this prefix. By default, `agouti` will be used.

## Example


## Break-and-Continue


## Graph Visualization


## Contributors

Simo Zhang
Luting Zhuo

## Support

Please report any issues or questions here on github or by email to simozhan@indiana.edu
