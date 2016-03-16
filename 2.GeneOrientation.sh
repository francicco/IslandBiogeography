#!/bin/bash

MITODIR=MetaGenomes
CodFrac=../Coding_len_frac.dat
rm Consistet.oriented.cntgs.dat
touch Consistet.oriented.cntgs.dat

cd $MITODIR

awk '{ if ($2 > 0.375) print $1}' $CodFrac | while read contig
do
	cd $contig
	MtDNA.GeneOrientation.check.py ../../collGeneOrientation.dat result $contig >> ../../Consistet.oriented.cntgs.dat
	cd ..
done



