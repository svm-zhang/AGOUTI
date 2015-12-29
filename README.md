# **AGOUTI**: Annotated Genome Optimization Using Transcriptome Information

Table of Contents
=================

  * [<strong>AGOUTI</strong>: Annotated Genome Optimization Using Transcriptome Information](#agouti-annotated-genome-optimization-using-transcriptome-information)
    * [Overview](#overview)
    * [Features](#features)
    * [Obtain AGOUTI](#obtain-agouti)
    * [Get Started](#get-started)
    * [Prepare Inputs](#prepare-inputs)
      * [Genome Assembly](#genome-assembly)
      * [SAM/BAM File](#sambam-file)
      * [Gene Models](#gene-models)
    * [Understand Outputs](#understand-outputs)
    * [Example](#example)
    * [Example Data](#example-data)
    * [Shred Assembly](#shred-assembly)
    * [Shred Annotation](#shred-annotation)
    * [Scaffoldding on Shredded Assembly](#scaffoldding-on-shredded-assembly)
    * [Break-and-Continue](#break-and-continue)
    * [Graph Visualization](#graph-visualization)
    * [Contributors](#contributors)
    * [Support](#support)

## Overview
AGOUTI uses paired-end RNA-seq reads to guide genome scaffolding and improve gene annotation. It works in the following steps:

1. Extracting uniquely mapped joining-read pairs
2. Denoise the set of joining-pairs using gene models
3. Traversing the graph built from the noise-free data to identify scaffolding paths
4. Updating assembly and gene annotation given the scaffolds

## Features

1. Scaffold hundreds to thousands of contigs, yielding more contiguous assemblies;
2. Reduce the number of gene models and update them simultaneously;
3. Shred genome assembly and corresponding annotation simultaneously;
4. Record any inconsistencies with the original (input) scaffolding results;
5. Recover original scaffolds from contigs not scaffolded by AGOUTI;
6. Support break-and-continue feature such that some time-consuming steps can be skipped if the previous run is successful;
7. Output graph in dot format for visualization

## Obtain AGOUTI

AGOUTI is designed to be light-weight by only requiring SAMtools, python 2.7 or above, and git 1.8.5 or later.

To download AGOUTI, please use:

    git clone https://github.com/svm-zhang/AGOUTI.git

To get the current version of the download, simply go to AGOUTI folder and do

    python agouti.py --version

You should see, for example, something like this:

    AGOUTI v0.2.4-dirty

Since v0.2.4, you can check if there is a new version by running:

    python agouti.py update

**It is RECOMMENDED to run update before you run AGOUTI.**

If you wish to use a specific version of AGOUTI, you can first fetch all the versions available:

    git fetch --all

Then show all available versions using:

    git tag

Next you simply need to specify a version. For example, if you'd like to use v0.2.3:

    git checkout v0.2.3 -b v0.2.3

You can also click on releases from AGOUTI github page to see all available versions.

## Get Started

In its simplest usage, AGOUTI takes three inputs: an initial genome assembly in FASTA format, paired-end RNA-seq reads mapped against the assembly in BAM format, and gene predictions from the initial assembly in GFF3 format. For instance:

    python agouti.py scaffold \
    -assembly example.fasta \
    -bam example.bam \
    -gff example.gff \
    -outdir ./example

This will produce a scaffoled assembly in FASTA format, and a updated gene models in GFF3 format. All files (including the intermediate files) will be stored under a directory specified by `-outdir`, "example" in this case.

## Prepare Inputs

### Genome Assembly

AGOUTI accepts assemblies as both contigs and scaffolds. In its scaffold form, AGOUTI breaks assemblies at gaps of a minimum lengths, essentially producing a shredded/split assembly (see **Shred Assembly**). AGOUTI scaffolds on the split assembly, and report any inconsistencies between the RNA-based scaffolding and the original scaffolding.

To shred a given assembly at gaps of at least 25 bp:

    python agouti.py shred \
    -assembly example.fasta \
    -p example \
    -mlg 25

This produces a shredded assembly: `example.ctg.fasta`, and a file of a format similar to Fasta: `example.shred.info.txt`. Each header gives IDs of sequences in the original assembly. Under each header is a list of pairs of shredded contigs and the length of gaps between them. A sequence without any gaps will be by itself, and NA are used for such cases.

**It is very important to use this split assembly in the following reads-mapping and gene prediction.**

### SAM/BAM File

Assuming you have a dataset of paired-end RNA-seq reads, `example.1.fq` and `example.2.fq`, and an assembly generated from either an assembler of your favorite or shredded by AGOUTI, `example.fasta` or `example.ctg.fasta`. You will first need to map the RNA-seq data against the assembly using a short-reads mapper, such as Bowtie2 or BWA. For example,

    bwa index example.fasta
    bwa mem -M example.fasta example.1.fq example.2.fq | samtools view - > example.bam

This produces a mapping results in BAM format. AGOUTI uses this BAM file for scaffolding. More specifically, it reads the file and extracts joining-pairs. A joining-pair is defined as one with both ends mapped to different contigs. AGOUTI uses only uniquely mapped ones by checking mapping quality. Short-reads mappers such as BWA, Bowtie2 uses a non-zero mapping quality to define unique mapping. If the mapper you are using does not use quality to mark ambiguous mapping, then you must first process your SAM/BAM file before running AGOUTI.

**Several more things worth of noting:**

1. Please run samtools flagstat to get stats of the mapping, and looks particular for number of pairs mapped to different chromosomes. If none, then AGOUTI will not be able to do any scaffolding.
2. Please make sure the BAM is sorted by reads name, not coordiantes.

### Gene Models

To run AGOUTI, you will also need a set of gene models predicted from the assembly. For instance,

    augustus --AUGUSTUS_CONFIG_PATH=[path to augustus config file] -gff3=on --species=[your sepcies] example.fasta > example.gff

At the end of gene prediction, you will now have a set of gene models predicted from the assembly. You can choose any * ab inito * gene predictor as long as it spits out the models in GFF format. More specifically, AGOUTI looks for the following information:

* lines annotated as `gene`
    * contig ID
    * gene ID, e.g. `ID=gene1` from the attribute column (i.e. last column)
    * start and stop positions of the gene
    * strand
* lines annotated as `CDS`
    * start and stop positions

**Several more things worth of noting:**

1. If your GFF file has FASTA sequences at the end (e.g. generated from MAKER pipeline), please make sure to use verions v0.2.5 or above.
2. It is also important to note that if there are no gene models found in your GFF, AGOUTI will issue an error and stop.

## Understand Outputs

AGOUTI outputs its results to a base directory specified by `-outdir`. Under the base director, there are several sub-folders created, each corresponding to a step built in AGOUTI. A run of AGOUTI using the command-line setting demonstrated in **Getting Started** will generated a structured output as shown in the following screenshot:

![Alt text](/image/example_outdir.png?raw=true "example output directory")

Each subfolder includes three types of file:

1. general progress meter info
2. debug info
3. intermediate outputs

To get a file with debug info you will need to specify `-debug`. An intermediate file can have all the joining-pairs, the denoised set of joining-pairs, the graph in DOT format, etc. Some intermediate files are important to support the break-and-continue feature, e.g. the file with the noise-free set of joining-pairs (see below for more details).

The `agouti.main.log` is prefixed with the string specified by `-p`, so do all the other files generated by AGOUTI. The sequence ID in the final assembly will also be as this prefix. By default, `agouti` will be used.

**The final assembly** and **the updated gene models** can be found under the base directory, `example`, along with plain text files of useful information, such as scaffolding paths, gene paths, differences between scaffolds generated by AGOUTI and original scaffolding.

## Example
Scaffoldding using joining-pairs with a minimum mapping quality of 20, a maximum of 5% mismtaches:

    python agouti.py scaffold -assembly example.fasta -bam example.bam -gff example.gff -outdir ./example -minMQ 20 -minFracMM 0.05

Scaffolding without updating gene model (v0.3.0 or above):

    python agouti.py scaffold -assembly example.fasta -bam example.bam -gff example.gff -outdir ./example -no_update_gff

Scaffolding a shredded assembly and report any inconsistencies between RNA-seq based scaffolding and orignial scaffolding:

    python agouti.py scaffold -assembly example.shred.fasta -bam example.bam -gff example.gff -outdir ./example -shredpath example.shred.info.txt

Shredding an assembly and annotation simultaneously (v0.3.0 or above):

    python agouti.py shred -assembly example.fasta -gff example.gff -p example

## Example Data

Here gives one [example data](http://www.indiana.edu/~hahnlab/software.html) set that we used in our paper.

## Shred Assembly

Given an assembly in its scaffold form, AGOUTI can shred scaffolds into contigs at gaps of a minimum length (5 by default, user-tunable). The following figure gives an example of how it works. Let's say a scaffold called `scaffold 1` in the assembly. This scaffold consists of three stretches of gaps of various lengths, 5, 3, and 9, respectively. By default, AGOUTI shreds it into three contigs, `Scaffold_1_0`, `Scaffold_1_1`, and `Scaffold_1_2`. AGOUTI does not cut at the second gap because it has a length of 3. Notably, AGOUTI uses `SEQID_INDEX` to tell the order of contigs in the given original scaffold. For scaffolds without gaps, AGOUTI does not split them.)

![Alt text](/image/shred_assembly.png?raw=true "example output directory")

## Shred Annotation

Since v0.3.0, AGOUTI is also able to shred gene annotation matching the give assembly. It compares start and end positions of features with coordinates of cut sites, and updates annotation accordingly. There are five types of features AGOUTI cares: gene, exon, CDS, firve_prime_UTR, three_prime_UTR. The following figure gives an example of how it works. Let's use the same scaffold (i.e. `Scaffold_1`) shredded in the picture above. Assume that there is a gene span across the second cut site, and it consists of three exons (green box). AGOUTI splits the assembly such that one gene becomes two (boxes in different colors) sitting on two different contigs. AGOUTI assigns them with different ID, in a similar fashion as names of shredded contigs, `GENEID_INDEX`. This naming tells 1) whether two shredded genes belong to a single one; and 2) the order.

![Alt text](/image/shred_annotation_1.png?raw=true "example output directory")

## Scaffoldding on Shredded Assembly

More details coming soon

## Break-and-Continue

AGOUTI is built with a couple of modules. The output of current module will be taken as the input as the next module. Given the same input, modules such as extracting joining-pairs from BAM file, spits out the same intermediate results. AGOUTI therefore tries to save some running time by skipping such steps if they were finished successfully from previous runs. To use this feature, simply run AGOUTI the second time with the same output directory and output prefix as the previous run. If you desire a fresh start, simply use `-overwrite` to overwrite all results generated previously, or gives a new prefix.

## Graph Visualization

AGOUTI makes the scaffolding graph accessible to users. Under `scaffolding` folder, you can find a file named after `[prefix].agouti_scaffolding.graph.dot`. The dot file can be directly loaded in packages like Graphviz. In the graph, contigs/vertices are in black circle, while there are two color codings for edges. Ones in red are the scaffolding path in the final assembly, and others in black are simply edges that were not traversed. Edges in dotted style represent connections with a number of supporting joining-pairs lower than the minimum specified.

## Contributors

Simo Zhang  
Luting Zhuo  
Matthew Hahn

## Support

Please report any issues or questions here on github or by email to simozhan@indiana.edu
