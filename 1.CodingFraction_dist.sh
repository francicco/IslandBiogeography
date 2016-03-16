#!/bin/bash

MITODIR=MetaGenomes
cd $MITODIR

rm tmp

ls -1 | while read contig
do
	cd $contig
	export SCF_len=`grep -v '^>' sequence.fas | awk '{sum+=(length$1)} END {print sum}'`
	export ANN_len=`grep gene result | cut -f 5,6 | awk '{sum+=($2-$1)} END {print sum+0.00000000001}'`
	#echo $contig, $SCF_len, $ANN_len
	echo $contig $ANN_len $SCF_len | awk '{frac=$2/$3} END {print $1, frac}' >> ../../tmp
	cd ..
done

cd ..

cat tmp | sort -k2,2rg > Coding_len_frac.dat

rm tmp


