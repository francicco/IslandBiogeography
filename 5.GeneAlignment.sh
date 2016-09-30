#!/bin/bash

WD=/scratch/c770/c7701100/Collembola_project/MITOS_server/tmp

allgeneFas=$WD/mtDNA.All.genes.fasta
ref_genes=$WD/RefCollMtDNA.formatted.fasta
MAC_DIR=/scratch/c770/c7701100/Collembola_project/MITOS_server
GAP_PERCENTAGE_CODING=33
GAP_PERCENTAGE_NONCOD=75
gene2sort='nad2 cox1 rrnL cox2 atp8 atp6 cox3 nad3 nad5 nad4 nad4l nad6 cob nad1 rrnL rrnS'
OVERWRITE=yes
KEEPTMP='no'

for gene in $gene2sort # Here the script process all genes one at the time
do 

	echo -e "\n#######################################################################"
	echo -e "\n`date` Align gene $gene."

	if [[ $OVERWRITE == 'yes' ]]
	then
		rm mtDNA.All.$gene.alned.NT.fasta
	fi

	if [ -f mtDNA.All.$gene.alned.NT.fasta ] # If there is this aligned fasta create it's consensus and deletes gaps otherwise it does the alignment.
	then

		echo -e "Generate consensus"
		AmbiguityConsensus.py -l $gene -i $WD/mtDNA.All.$gene.alned.NT.fasta -o $WD/mtDNA.All.$gene.consensus.fasta
		
		echo -e "Remove gaps from consensus"
		GapAlnRemover.py -c $WD/mtDNA.All.$gene.consensus.fasta -i $WD/mtDNA.All.$gene.alned.NT.fasta -t coding -o $WD/mtDNA.All.$gene.alned.NT.ungapped.fasta
	else
		ls $allgeneFas 
		# Extract single gene from file containint all genes and create file to align. Creating the dataset to align.
		grep --no-group-separator -A1 -w $gene $allgeneFas > $WD/mtDNA.All.$gene.fasta
		grep --no-group-separator -A1 -w $gene $ref_genes >> $WD/mtDNA.All.$gene.fasta
	
		# Count file length
		Aln_len=`grep -c '^>' $WD/mtDNA.All.$gene.fasta`

		# Start iteration as long as the length of the file.
		for i in $(seq $Aln_len)
		do
			echo -e "\n`date` Iteration $i started."
		
			iter=`bc <<< $i+1`

			if [ $i == 1 ]
			then
				ALNMENT=$WD/mtDNA.All.$gene
			else
				ALNMENT=$WD/mtDNA.All.$gene.$i.Iteration
			fi

			N_SEQ=`grep -c '^>' $ALNMENT.fasta`
			echo -e "`date` Alignment of $N_SEQ sequences in $ALNMENT.fasta file."

			if [[ $gene == rrnL || $gene == rrnS ]]
			then
				GEPTR=$GAP_PERCENTAGE_NONCOD
				sed -ibk 's/;//g' $ALNMENT.fasta
				rm $ALNMENT.fastabk

				cat $ALNMENT.fasta

				clustalw2 -infile=$ALNMENT.fasta -outfile=$ALNMENT.alned.fasta \
					-output=FASTA -dnamatrix=IUB \
					-gapopen=10 -gapext=0.2 \
					-gapdist=10 -iteration=ALIGNMENT \
					-numiter=1000 -clustering=UPGMA
			else
				GEPTR=$GAP_PERCENTAGE_CODING
				java -jar $WD/macse_v1.01b.jar \
        			     -prog alignSequences \
        			     -gc_def 5 \
        			     -out_NT $ALNMENT.alned.fasta \
        			     -seq $ALNMENT.fasta
			fi

			cat $ALNMENT.alned.fasta
	
			echo -e "\n`date` Check and exclude short sequences from the alignment." 
			echo -e "gapCount.py $ALNMENT.alned.fasta $GEPTR tmp | cut -d ' ' -f 1 | sort | cut -f 1| sed 's/_H_/; /g' | cut -d ' ' -f 1 > $WD/mtDNA.All.$gene.alned.$i.ids"
			gapCount.py $ALNMENT.alned.fasta $GEPTR $WD/mtDNA.All.$gene.$iter.Iteration.fasta > $WD/mtDNA.All.$gene.alned.$i.ids
			
			p_iter=`bc <<< $i-1`
			
			if [ -f $WD/mtDNA.All.$gene.alned.$p_iter.ids ]
			then
				FILE1=`wc -l $WD/mtDNA.All.$gene.alned.$p_iter.ids | cut -d ' ' -f 1`
				FILE2=`wc -l $WD/mtDNA.All.$gene.alned.$i.ids | cut -d ' ' -f 1`
				if [[ $FILE1 == $FILE2 || $FILE2 == 0 ]]
				then
					if [[ $gene == rrnL || $gene == rrnS ]]
					then
						echo -e "Removing temporary files"
                                        	sed 's/_H_/; /g' $ALNMENT.alned.fasta > $WD/mtDNA.All.$gene.alned.NT.fasta
					else
						echo -e "Removing temporary files"
						mv $ALNMENT.alned.fasta $WD/mtDNA.All.$gene.alned.NT.fasta
						mv $WD/mtDNA.All.$gene.$i.Iteration_macse_AA.fasta $WD/mtDNA.All.$gene.alned.AA.fasta
					fi

					echo -e "Generate consensus"
					AmbiguityConsensus.py -l $gene -i $WD/mtDNA.All.$gene.alned.NT.fasta -o $WD/mtDNA.All.$gene.consensus.fasta
					echo -e "Remove gaps from consensus"
					GapAlnRemover.py -c $WD/mtDNA.All.$gene.consensus.fasta -i $WD/mtDNA.All.$gene.alned.NT.fasta -t coding -o $WD/mtDNA.All.$gene.alned.NT.ungapped.fasta
					
					if [[ $KEEPTMP='no' ]]
					then
						rm -f $WD/tmp $WD/*.ids $WD/*.Iteration.* $WD/*.alned.fasta $WD/*_macse_AA.fasta $WD/mtDNA.All.$gene.fasta $WD/*.dnd
					fi
					break
				fi
			fi
	
			n_long=`wc -l $WD/mtDNA.All.$gene.alned.$i.ids | cut -d ' ' -f 1`
			n_short=`bc <<< $Aln_len-$n_long`
		
		done

	fi
done


echo -e "FASTA alignment conversion to NEXUS and concatenation"
ConcatenateFasta2Nexus.py -i $WD -s .ungapped.fasta -o $WD/mtDNA.All.concatenated.nex

