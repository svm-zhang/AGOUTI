# **AGOUTI**: Annotated Genome Optimization Using Transcriptome Information

## Overview
AGOUTI uses paired-end RNA-seq reads to guide genome scaffolding and improve gene annotation. It works in the following steps:

1. Extracting uniquely mapped joining-read pairs
2. Denoise the set of joining-pairs using gene models
3. Traversing the graph built from the noise-free data to identify scaffolding paths
4. Updating assembly and gene annotation given the scaffolds

## Features

1. Scaffold hundreds to thousands of contigs, yielding more contiguous assemblies;
2. Reduce the number of gene models and update them simultaneously;
3. Record any inconsistencies with the original (input) scaffolding results;
4. Support break-and-continue feature such that some time-consuming steps can be skipped if the previous run is successful;
5. Output graph in dot format for visualization

## Obtain AGOUTI

To download AGOUTI, please use git to download the most recent development version through:

    git clone https://github.com/svm-zhang/AGOUTI.git

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
usage: agouti.py [-h] ...

Welcome to AGOUTI!

optional arguments:
    -h, --help  show this help message and exit

Commands:

    scaffold    Scaffolding genome assembly and update genome annotation
    shred       Shredding genome assembly into contigs at gaps of a minimum length
```
### Scaffold interface

To run scaffolding using AGOUTI, simply:

    python agouti.py scaffold -h

```
usage: agouti.py [-h] -contig FILE -bam FILE -gff FILE -out DIR [-mnl INT] [-p STR]

Welcome to AGOUTI!

optional arguments:
	-h, --help                  show this help message and exit
	-assembly           FILE    specify the assembly in FASTA format
	-bam                FILE    specify the RNA-seq mapping results in BAM format
	-gff                FILE    specify the predicted gene model in GFF format
	-outdir             DIR     specify the base directory to store all output files
	-p                  STR     specify the prefix for all output files [agouti]
	-k                  INT     minimum number of joining reads pairs supports [5]
    -nN                 INT     number of Ns put in between a pair of contigs [1000]
    -minMQ              INT     minimum mapping quality to use [5]
    -minFracOvl         FLOAT   minimum percentage of alignment length: alnLen/readLen [0.0]
    -maxFracMM          FLOAT   maximum fraction of mismatch per alignment [1.0]
    -debug                      specify to have info for debugging
    -overwrite                  specify to overwirte all results from previous run
    -v, --version               show program's version number and exit
```

## Get Started

In it simplest usage, AGOUTI takes three inputs: an initial genome assembly in FASTA format, paired-end RNA-seq reads mapped against the assembly in BAM format, and gene predictions from the initial assembly in GFF3 format. For instance:

    python agouti.py \
    -assembly example.fasta \
    -bam example.bam \
    -gff example.gff \
    -outdir ./example

This will produce a scaffoled assembly in FASTA format, and a updated gene models in GFF3 format. All files (including the intermediate files) will be stored under a directory specified by `-outdir`, "example" in this case.

## Prepare Inputs

Assuming you have a dataset of paired-end RNA-seq reads, `example.1.fq` and `example.2.fq`, and an initial assembly generated from an assembler of your favorite, `example.fasta`. You will first need to map the RNA-seq data against the assembly using a short-reads mapper, such as BWA. For example,

    bwa index example.fasta
    bwa mem -M example.fasta example.1.fq example.2.fq | samtools view - > example.bam

At the end of reads mapping, you will have the mapping results in BAM format. AGOUTI uses only uniquely mapped joining-pairs by checking mapping quality. For short-reads mappers such as BWA, you can use mapping quality to tell unique or not, e.g. mapQ > 0. If you use other mappers, you can filter out the ambiguous ones prior to input to AGOUTI. In the future version, we will make AGOUTI being able too internally recognize ambiguous reads mapping from different mappers.

AGOUTI expects each reads pair to be next to each other in the SAM/BAM file. Therefore, there is no need to sort the BAM file by coordiantes. In addition, please FILTER secondary and supplementary alignments from the SAM/BAM file. SAMTOOLs can do the filter by specifying -F options.

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

To get a file with debug info you will need to specify `-debug`. An intermediate file can have all the joining-pairs, the denoised set of joining-pairs, the graph in DOT format, etc. Some intermediate files are important to support the break-and-continue feature, e.g. the file with the noise-free set of joining-pairs (see below for more details).

The `agouti.main.log` is prefixed with the string specified by `-p`, so do all the other files generated by AGOUTI. The sequence ID in the final assembly will also be as this prefix. By default, `agouti` will be used.

**The final assembly** and **the updated gene models** can be found under `agouti_update`, along with files of useful information, such as the mapping between new scaffold and old contigs.

## Break-and-Continue

AGOUTI is built with a couple of modules. The output of current module will be taken as the input as the next module. Given the same input, modules such as extracting joining-pairs from BAM file, spits out the same intermediate results. AGOUTI therefore tries to save some running time by skipping some steps if they were finished successfully from previous runs. To use this feature, simply run AGOUTI the second time with the same output prefix as the previous run. In some cases you want to have a fresh start, simply use `-overwrite` to overwrite all results generated previously, or gives a new prefix.

## Graph Visualization

AGOUTI makes the scaffolding graph accessible to users. Under `scaffolding` folder, you can find a file named after `[prefix].agouti_scaffolding.graph.dot`. The dot file can be directly loaded in packages like Graphviz. In the graph, contigs/vertices are in black circle, while there are two color codings for edges. Ones in red are the scaffolding path in the final assembly, and others in black are simply edges that were not traversed. Edges in dotted style represent connections with a number of supporting joining-pairs lower than the minimum specified.

## Contributors

Simo Zhang  
Luting Zhuo  
Matthew Hahn

## Support

Please report any issues or questions here on github or by email to simozhan@indiana.edu
