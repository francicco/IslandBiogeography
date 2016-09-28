#!/bin/bash


MITODIR=/scratch/c7701100/Collembola_project/MITOS_server/MetaGenomes
Good_scf=/scratch/c7701100/Collembola_project/MITOS_server/results_100contigs

touch annotations.bed

cd $MITODIR

ls -1 | while read contig
do
	cd $contig
	cut -f 1-8 result >> ../../annotations.bed
	cd ..
done


cd $MITODIR

touch annotations.bed.good

cut -d ':' -f1 $Good_scf | while read contig
do
        cd $contig
        cut -f 1-8 result >> ../../annotations.bed.good
        cd ..
done
