#!/bin/bash

WD=/scratch/c7701100/Collembola_project/MITOS_server

MITODIR=MetaGenomes
mtDNAscf=$WD/mtDNA.All.scaffolds.ids
allscfFas=$WD/mtDNA.All.scaffolds.fasta
allscfBed=$WD/mtDNA.All.scaffolds.bed
allgeneFas=$WD/mtDNA.All.genes.fasta


rm $allscfFas $allscfBed $allgeneFas

touch $allscfFas
touch $allscfBed
touch $allgeneFas

cd $MITODIR

cut -d ':' -f 1 $mtDNAscf | while read scf
do
	cd $scf
	echo Proessing $scf
	echo Conversion...
	MtDNA.mitos_result2bedfile.py result > $scf.ann.bed
	echo Extract fasta seq...
	MtDNA.Bed2fasta.py -f sequence.fas -b $scf.ann.bed -o $scf.genes.fasta
	echo -e "Concatenate results...\n"
	cat sequence.fas >> $allscfFas
	cat $scf.ann.bed >> $allscfBed
	cat $scf.genes.fasta >> $allgeneFas
	cd ..
done




