#!/bin/bash

WD=/scratch/c7701100/Collembola_project/MITOS_server

MITODIR=MetaGenomes
Cntgs_oriented=$WD/Consistet.oriented.cntgs.dat
GeneOrder=$WD/collGeneOrder.dat

tmp=$WD/.tmp
single=$WD/results_snglGene.contigs
All=$WD/results_All.contigs
Allmin1=$WD/results_All-1.contigs

rm $tmp
touch $tmp

cd $MITODIR

cut -d ' ' -f 1 $Cntgs_oriented | while read scf
do
	cd $scf
	export strand=`grep $scf $Cntgs_oriented | cut -d ' ' -f 2`
	#echo $scf $strand
	MtDNA.GeneOrder.check.py cntg_g_order.dat $GeneOrder $scf $strand >> $tmp
	MtDNA.GeneOrder.check.py cntg_g_order.dat $GeneOrder $scf $strand | grep 'only 1 gene' >> $single
	MtDNA.GeneOrder.check.py cntg_g_order.dat $GeneOrder $scf $strand | awk '{ if ($2 == $4) print $0}' >> $All
	MtDNA.GeneOrder.check.py cntg_g_order.dat $GeneOrder $scf $strand | awk '{ if ($2-2 == $4) print $0}' >> $Allmin1
	cd ..
done

cat $tmp | grep 'only 1 gene' > $single
cat $tmp | awk '{ if ($2 == $4) print $0}' > $All
cat $tmp | awk '{ if ($2-2 == $4) print $0}' > $Allmin1



