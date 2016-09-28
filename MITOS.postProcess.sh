#!/bin/bash


MITODIR=MetaGenomes
cd $MITODIR

#ls -1 | while read contig
#do
#	unzip $contig
#	rm $contig
#done


ls -1 | while read directory
do
	cd $directory
	touch $directory
	contig=`sed -n 1p sequence.fas | sed 's/^>//'`
	cd ..
	mv $directory $contig
done


