#!/bin/bash


#####################################
#### Filter putative mtDNA contigs obtained in step 1 (1.CodingFraction_dist) based on gene orientation. 
#### Takes as input a directory containing all the files generated by the MITOS webserver and in step 1.    
####
#### Usage: bash 2.GeneOrientation.sh Directory Prefix Output  
####        - Directory = directory with all MITOS raw data folders                              
####        - Prefix = ID of the sample being processed
####        - Output = Path to the output folder of your choice (should be the same as in step 1)
####
#### Dependencies: MtDNA.GeneOrientation.check.py (in bin folder), collGeneOrientation.dat file (in bin folder).
####
#### The collGeneOrientation.dat is a space delimited file containing the mitochondrial genes (coding and RNA)
#### in the good order for the species/genus of interest (check in METAMIGA for example) but with the formatting as in 
#### MITOS (for example: ND2 would be nad2 to fit the MITOS annotations). Each gene name is followed by "+" or "-" 
#### according to whether the gene is on the positive of negative strand. 
#### E.G: nad2+ cox1+ cox2+ etc... 
####
#### Francesco Cicconardi (2016)
#####################################


### Getting script arguments
MITODIR=$1
PREFIX=$2
OUTPUT=$3                    

### Important parameter
CODFRAC=$OUTPUT/$PREFIX.Coding_len_frac.dat    # File generated in step 1

### Removing old files
rm $OUTPUT/$PREFIX.Consistet.oriented.cntgs.dat
touch $OUTPUT/$PREFIX.Consistet.oriented.cntgs.dat

### Looping over all contigs
cd $MITODIR

awk '{ if ($2 > 0.375) print $1}' $CODFRAC | while read contig
do
    cd $contig
    python ../../bin/MtDNA.GeneOrientation.check.py ../../bin/collGeneOrientation.dat result $contig >> $OUTPUT/$PREFIX.Consistet.oriented.cntgs.dat
    cd ..
done


