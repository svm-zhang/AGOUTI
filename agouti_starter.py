#!/usr/bin/python

"""
Create a bash file for running AGOUTI pipeline

last modified: 10-9-2014

"""

import os
import sys
import re
import argparse
import subprocess 

def main():
    
    args = run_argparser()
   
    contigFA = args.contigFA
    inBam = args.inBam
    outDir = args.outDir
    sampleID = args.sampleID
    sspacePATH = args.sspacePATH
    agoutiPATH = args.agoutiPATH
    isPE = not args.isSE    
    
    print outDir
    # to-do: load the configuration file, such as settings for sspace program

    
    # create tmp folder under outDir
    moveOn = False
    if os.path.isdir(outDir):
        moveOn = True
        tmpDir = os.path.join(outDir,"tmp")
        if not os.path.isdir(tmpDir):
            os.mkdir(tmpDir)
    if not  moveOn:
        sys.exit("please double check if sspacePATH is valid")
    
    # get inter-contig discordant read pairs or split-reads from input Bam
    if isPE:
       getDiscdtReads = get_discdt_pe(inBam,tmpDir,sampleID,agoutiPATH) 
    
    else :
       getDiscdtReads = get_discdt_se(inBam,tmpDir,sampleID,agoutiPATH)
    
    # perfrom RNAseq-guided scaffolding using SSPACE(STANDARD version)
    runSSPACE = run_sspace(contigFA,agoutiPATH,sspacePATH)
    
    # now parse SSPACE output since this program is desgined for DNAseq-guided scaffolding
    parseSSPACE = parse_sspace(outDir)
    
    # to-do: create an output .bash file, might add PBS header if necessary
    
    
    


# end of main


def get_discdt_pe(inBam,tmpDir,prefix,agoutiPATH):
    """
    create commands that extract discordant read pairs from input bam file,
    reads with mapping quality <5 will be removed.
    """
    command = """
prefix=%s
inBam=%s
agoutiPATH=%s
tmp=%s

echo start extracting discordantly mapped read pairs from input bam file...

# extract inter-contig discordant read pairs, make sure that bam file has been sorted by read name
samtools view -F 268 -q 5  $inBam |awk '{ if($7!="*" && $7!="=") print $0}'> $tmp/$prefix.q5.aln


# remove singletons
python $agoutiPATH/getDiscdtPE.py  $tmp/$prefix.q5.aln  $tmp/$prefix.discdt.q5.aln

# convert back to bam format for next step
samtools view -H $inBam  |cat - $tmp/$prefix.discdt.q5.aln | samtools view -bS   -  > $tmp/$prefix.discdt.q5.bam

# create .bedpe file
bamToBed  -bedpe  -i $tmp/$prefix.discdt.q5.bam  > $tmp/$prefix.discdt.q5.bedpe

# create .linkpair (linked-contig) file for SSPACE
cat $tmp/$prefix.discdt.q5.bedpe | python $agoutiPATH/bedpe2sspacelinkpair.py > $tmp/$prefix.sspace.linkpair

"""%(prefix,inBam,agoutiPATH,tmpDir)
    
    print command
    return command


def get_discdt_se(inBam,tmpDir,prefix,agoutiPATH):
    """
     create commands that extract discordant read pairs from input bam file,
     reads with mapping quality <5 will be removed.
    """
    command = """
prefix=%s
inBam=%s
agoutiPATH=%s
tmp=%s

echo start extracting split-read from input bam file...

# extract inter-contig discordant split reads, make sure that bam file has been sorted by read name
samtools view -q 5 $inBam |awk '{if(match($0,"SA:Z")) print $0}' > $tmp/$prefix.split.aln
python $agoutiPATH/getDiscdtSE.py $tmp/$prefix.split.aln >  $tmp/$prefix.discdt.q5.txt

# convert back to bam format for next step
samtools view -H $inBam |cat - $tmp/$prefix.discdt.q5.txt > $tmp/$prefix.discdt.q5.sam

# create a linkpair (linked-contig) file in .bedpe format
python $agoutiPATH/makeBedpeFromSE.py $tmp/$prefix.discdt.q5.sam  $tmp/$prefix.discdt.q5.sam.se.bedpe    

# create .linkpair (linked-contig) file as the input for SSPACE 
cat $tmp/$prefix.discdt.q5.sam.se.bedpe | python $agoutiPATH/bedpe2sspacelinkpair.py > $tmp/$prefix.sspace.linkpair

"""%(prefix,inBam,agoutiPATH,tmpDir)

    print command
    return command

def run_sspace(contigFA,agoutiPATH,sspacePATH):
    """
    run SSPACE-STANDARD to perform scaffolding
    """
    command = """

contigFA=%s
sspacePATH=%s

echo start running SSPACE ...

# create Lib file storing the library for SSPACE scaffolding
# Noted that in Lib_ file, the  mean distance (4th column) between the
# paired reads (or split-read) and standard deviation (5th column)
# is set  to a much higher value than the values used by DNA seq-guided scaffolding.
# Such setting is to facilitate the intron-depleted nature of RNAseq reads.
 
echo "Lib1 TAB  $tmp/$prefix.sspace.linkpair  100000   0.9  FR" > $tmp/Lib_$prefix

cd $tmp
# run SSPACE STANDARD version
perl  $sspacePATH/SSPACE_Standard_v3.0.pl \
-l $tmp/Lib_$prefix \
-s $contigFA \
-k 2 -x 0 -z 0 -a 0.7 -n 10 -T 8 -p 0 \
-b $prefix.sspaceout

cd -
"""%(contigFA,sspacePATH)
    
    print command
    return command


def parse_sspace(outDir):
    """
    Parse the output of SSPACE into scaffolding paths:
    1) convert mis-calculated gap into fix-sized gap (100'N')
    2) generate scaffolding paths
    """
    command = """
sspaceoutDir=$tmp/$prefix.sspaceout
outDir=%s

# parse the output of sspace to convert incorrectly calculated gap size (since SSPACE assume input to be mate-pair lib) into a string of 100 "N"s
python $agoutiPATH/convertgap_SSPACEOUT.py  $sspaceoutDir/$prefix.sspaceout.final.scaffolds.fasta $tmp

# add scaffolding paths to final output
grep ">" $sspaceoutDir/intermediate_results/$prefix.sspaceout.formattedcontigs*fasta > $tmp/$prefix.sspaceout.keys

python  $agoutiPATH/parseSSPACEOUT--forBenchMark.py \
    $tmp/gapparsed-$prefix.sspaceout.final.scaffolds.fasta   \
    $sspaceoutDir/$prefix.sspaceout.final.evidence  \
    $tmp/$prefix.sspaceout.keys  \
    $outDir/$prefix.agouti.fasta

# compute N50
python $agoutiPATH/calculateN50.py  $outDir/$prefix.agouti.fasta
"""%(outDir)

    print command
    return command




def run_argparser():
    """
    """
    
    use_message = '''
    
    Welcome to AGOUTI!\n
    Before running this script, make sure you have these software installed:\n
        SAMtools\n
        BEDtools\n
        SSPACE(STANDARD version)\n

    '''
    
    parser = argparse.ArgumentParser(description=use_message)

    parser.add_argument("-contig",
                        dest="contigFA",
                        required=True,
                        help="specify the contig fasta file to be scaffolded")
    parser.add_argument("-bam",
                        dest="inBam",
                        required=True,
                        help="specify the .bam file generated by BWA mem")
    parser.add_argument("-out",
                        dest="outDir",
                        default=".",
                        required=True,
                        help="specify the directory to store output files")
    parser.add_argument("-sampleID",
                        dest="sampleID",
                        default="agouti",
                        help="specify the sampleID as the prefix of output files")
    parser.add_argument("-isSE",
                        dest="isSE",
                        action="store_true",
                        default=False,
                        help="turn on this option if input is single-end RNAseq reads")
    parser.add_argument("-sspacePATH",
                        dest="sspacePATH",
                        default="../SSPACE-STANDARD-3.0_linux-x86_64",
                        help="specify the path to SSPACE")
    parser.add_argument("-agoutiPATH",
                        dest="agoutiPATH",
                        default=".",
                        help="specify the path to AGOUTI")
    # other possible arguments include "-compressFQ","-doFASTQC","doRSeqQC"
 
    args = parser.parse_args()
    return args





if __name__ == "__main__":
    main()
