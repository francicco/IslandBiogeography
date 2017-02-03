#!/bin/bash


#####################################
#### Extracts single fasta files from single genes.
#### Uses as input contigs from the previous step (3.GeneOrder).     
####
#### Usage: bash 4.GetAnnFasta.sh Directory Prefix Output
####        - Directory = directory with all MITOS raw data folders                               
####        - Prefix = ID of the sample being processed
####        - Output = Path to the output folder of your choice (must be the same as in steps 1,2 and 3)
####
#### Dependencies: MtDNA.mitos_result2bedfile.py (in bin folder), MtDNA.Bed2fasta.py (in bin folder).
####
#### The script generates three outputs: 
####    1) A fasta file containing the complete sequences of all the annotated contigs (All.contigs.fasta)
####    2) A bed file summarizing the position of each genes (All.contigs.bed)
####    3) A fasta file containing the sequence of each gene found (All.genes.fasta)
####
#### Francesco Cicconardi (2016)
#####################################

### Getting script arguments
MITODIR=$1                
PREFIX=$2
OUTPUT=$3 

### Important parameter
CTG_ID=$OUTPUT/$PREFIX.All.contigs.ids

### Defining output files to be generated
CTG_FAS=$OUTPUT/$PREFIX.mtDNA.All.contigs.fasta
CTG_BED=$OUTPUT/$PREFIX.mtDNA.All.contigs.bed
GENE_FAS=$OUTPUT/$PREFIX.mtDNA.All.genes.fasta

### Removing old files
rm $CTG_FAS $CTG_BED $GENE_FAS

touch $CTG_FAS
touch $CTG_BED
touch $GENE_FAS

### Entering directory
cd $MITODIR

### Keeping only the ID of the contigs and looping over them
cut -d ':' -f 1 $CTG_ID | while read contig
do
    cd $contig
    echo Processing $contig
    echo Conversion...
    python ../../bin/MtDNA.mitos_result2bedfile.py result > $contig.ann.bed
    echo Extract fasta seq...
    python ../../bin/MtDNA.Bed2fasta.py -f sequence.fas -b $contig.ann.bed -o $contig.genes.fasta
    echo -e "Concatenate results...\n"
    cat sequence.fas >> $CTG_FAS
    cat $contig.ann.bed >> $CTG_BED
    cat $contig.genes.fasta >> $GENE_FAS
    cd ..
done
