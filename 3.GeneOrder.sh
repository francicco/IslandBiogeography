#!/bin/bash


#####################################
#### Filter putative mtDNA contigs based on gene order.    
#### Uses as input contigs from the previous step (2.GeneOrientation). 
####
#### Usage: bash 3.GeneOrder.sh Directory Prefix Output
####        - Directory = directory with all MITOS raw data folders                               
####        - Prefix = ID of the sample being processed
####        - Output = Path to the output folder of your choice (must be the same as in steps 1 and 2)
####
#### Dependencies: MtDNA.GeneOrder.check.py (in bin folder), collGeneOrder.dat file (in bin folder).
####
#### The collGeneOrder.dat is a space delimited file containing the mitochondrial genes (coding and RNA)
#### in the good order for the species/genus of interest (check in METAMIGA for example) but with the formatting as in 
#### MITOS (for example: ND2 would be nad2 to fit the MITOS annotations). Genes on the positive strand are 
#### followed by "+". 
#### E.G: nad2+ cox1+ cox2+ nad1 etc... 
####
#### Francesco Cicconardi (2016)
#####################################

### Getting script arguments
MITODIR=$1                
PREFIX=$2
OUTPUT=$3 

### Defining output files to be generated
TMP=$OUTPUT/tmp.gene.order
SINGLE=$OUTPUT/$PREFIX.results_snglGene.contigs
ALL=$OUTPUT/$PREFIX.results_All.contigs
ALLMIN1=$OUTPUT/$PREFIX.results_All-1.contigs

### Removing old temporary files
rm $TMP
touch $TMP

### Changing directory
cd $MITODIR

### Keeping only the column with the name of the contigs from the previous step, and looping over them
cut -d ' ' -f 1 $OUTPUT/$PREFIX.Consistet.oriented.cntgs.dat | while read contig
do
    cd $contig
    export strand=`grep $contig $OUTPUT/$PREFIX.Consistet.oriented.cntgs.dat | cut -d ' ' -f 2`
    echo $contig $strand
    python ../../bin/MtDNA.GeneOrder.check.py cntg_g_order.dat ../../bin/collGeneOrder.dat $contig $strand >> $TMP
    python ../../bin/MtDNA.GeneOrder.check.py cntg_g_order.dat ../../bin/collGeneOrder.dat $contig $strand | grep 'only 1 gene' >> $SINGLE
    python ../../bin/MtDNA.GeneOrder.check.py cntg_g_order.dat ../../bin/collGeneOrder.dat $contig $strand | awk '{ if ($2 == $4) print $0}' >> $ALL
    python ../../bin/MtDNA.GeneOrder.check.py cntg_g_order.dat ../../bin/collGeneOrder.dat $contig $strand | awk '{ if ($2-2 == $4) print $0}' >> $ALLMIN1
    cd ..
done

### Generating the final output 
cat $SINGLE $ALL $ALLMIN1 > $OUTPUT/$PREFIX.All.contigs.ids
