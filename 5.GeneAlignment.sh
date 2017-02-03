#!/bin/bash


#####################################
#### Iterative alignment of protein-coding sequences with MACSE and of non-coding sequences with CLUSTAL-W
#### Uses as input gene sequences from the previous step (4.GetAnnFasta).      
####
#### Usage: bash 5.GeneAlignment.sh Prefix Output bin
####        - Prefix = ID of the sample being processed
####        - Output = Path to the output folder of your choice (must be the same as in steps 1,2,3 and 4)
####        - bin = Path for the bin folder
####
#### Dependencies: AmbiguityConsensus.py, GapAlnRemover.py, RefCollMtDNA.formatted.fasta, gapCount.py, ConcatenateFasta2Nexus.py, macse_v1.01b.jar (all in bin folder) and ClustalW2
####
#### The RefCollMtDNA.formatted.fasta file is a file containing all the gene sequences available of the
#### species/genus of interest dowloaded from a database (METAMIGA, GENBANK). 
####
#### Francesco Cicconardi (2016)
#####################################



### Getting script arguments
PREFIX=$1
OUTPUT=$2
BIN=$3

### Important parameters
allgeneFas=$OUTPUT/$PREFIX.mtDNA.All.genes.fasta #Output from the step 4
ref_genes=bin/RefCollMtDNA.formatted.fasta #Gene sequences from public database
GAP_PERCENTAGE_CODING=33 #Removing sequences shorter than the 33% of the overall alignment, in coding sequences
GAP_PERCENTAGE_NONCOD=75 #Removing sequences shorter than the 75% of the overall alignment, in non-coding sequences
gene2sort='nad2 cox1 rrnL cox2 atp8 atp6 cox3 nad3 nad5 nad4 nad4l nad6 cob nad1 rrnL rrnS'
OVERWRITE=yes     #If 'yes', overwrite previous run
KEEPTMP='no'     #If 'no', remove temporary files


### Alignement of coding sequences, all the genes at one time
for gene in $gene2sort 
do 

    echo -e "\n#######################################################################"
    echo -e "\n`date` Align gene $gene."

    if [[ $OVERWRITE == 'yes' ]]
    then
        rm $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta
    fi

    if [ -f $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta ] # If there is this aligned fasta create it's consensus and deletes gaps otherwise it does the alignment.
    then

        echo -e "Generate consensus"
        python bin/AmbiguityConsensus.py -l $gene -i $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta -o $OUTPUT/$PREFIX.mtDNA.All.$gene.consensus.fasta
        
        echo -e "Remove gaps from consensus"
        python bin/GapAlnRemover.py -c $OUTPUT/$PREFIX.mtDNA.All.$gene.consensus.fasta -i $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta -t coding -o $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.ungapped.fasta
    else
        ls $allgeneFas 
        # Extract single gene from file containint all genes and create file to align. Creating the dataset to align.
        grep --no-group-separator -A1 -w $gene $allgeneFas > $OUTPUT/$PREFIX.mtDNA.All.$gene.fasta
        #grep --no-group-separator -A1 -w $gene $ref_genes >> $OUTPUT/$PREFIX.mtDNA.All.$gene.fasta #Keeps reference genes in the final alignment
    
        # Count file length
        Aln_len=`grep -c '^>' $OUTPUT/$PREFIX.mtDNA.All.$gene.fasta`

        # Start iteration as long as the length of the file.
        for i in $(seq $Aln_len)
        do
            echo -e "\n`date` Iteration $i started."
        
            iter=`bc <<< $i+1`

            if [ $i == 1 ]
            then
                ALNMENT=$OUTPUT/$PREFIX.mtDNA.All.$gene
            else
                ALNMENT=$OUTPUT/$PREFIX.mtDNA.All.$gene.$i.Iteration
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
                java -jar bin/macse_v1.01b.jar \
                         -prog alignSequences \
                         -gc_def 5 \
                         -out_NT $ALNMENT.alned.fasta \
                         -seq $ALNMENT.fasta
            fi

            cat $ALNMENT.alned.fasta
    
            echo -e "\n`date` Check and exclude short sequences from the alignment." 
            echo -e "gapCount.py $ALNMENT.alned.fasta $GEPTR tmp | cut -d ' ' -f 1 | sort | cut -f 1| sed 's/_H_/; /g' | cut -d ' ' -f 1 > $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.$i.ids"
            python bin/gapCount.py $ALNMENT.alned.fasta $GEPTR $OUTPUT/$PREFIX.mtDNA.All.$gene.$iter.Iteration.fasta > $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.$i.ids
            
            p_iter=`bc <<< $i-1`
            
            if [ -f $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.$p_iter.ids ]
            then
                FILE1=`wc -l $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.$p_iter.ids | cut -d ' ' -f 1`
                FILE2=`wc -l $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.$i.ids | cut -d ' ' -f 1`
                if [[ $FILE1 == $FILE2 || $FILE2 == 0 ]]
                then
                    if [[ $gene == rrnL || $gene == rrnS ]]
                    then
                        echo -e "Removing temporary files"
                                            sed 's/_H_/; /g' $ALNMENT.alned.fasta > $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta
                    else
                        echo -e "Removing temporary files"
                        mv $ALNMENT.alned.fasta $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta
                        mv $OUTPUT/$PREFIX.mtDNA.All.$gene.$i.Iteration_macse_AA.fasta $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.AA.fasta
                    fi

                    echo -e "Generate consensus"
                    python bin/AmbiguityConsensus.py -l $gene -i $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta -o $OUTPUT/$PREFIX.mtDNA.All.$gene.consensus.fasta
                    echo -e "Remove gaps from consensus"
                    python bin/GapAlnRemover.py -c $OUTPUT/$PREFIX.mtDNA.All.$gene.consensus.fasta -i $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.fasta -t coding -o $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.NT.ungapped.fasta
                    
                    if [[ $KEEPTMP='no' ]]
                    then
                        rm -f $OUTPUT/$PREFIX.tmp $OUTPUT/$PREFIX.*.ids $OUTPUT/$PREFIX.*.Iteration.* $OUTPUT/$PREFIX.*.alned.fasta $OUTPUT/$PREFIX.*_macse_AA.fasta $OUTPUT/$PREFIX.mtDNA.All.$gene.fasta $OUTPUT/$PREFIX.*.dnd
                    fi
                    break
                fi
            fi
    
            n_long=`wc -l $OUTPUT/$PREFIX.mtDNA.All.$gene.alned.$i.ids | cut -d ' ' -f 1`
            n_short=`bc <<< $Aln_len-$n_long`
        
        done

    fi
done

cd $OUTPUT

echo -e "FASTA alignment conversion to NEXUS and concatenation"
python $BIN/ConcatenateFasta2Nexus.py -i $OUTPUT -s .ungapped.fasta -o $OUTPUT/$PREFIX.mtDNA.All.concatenated.nex

cd .. 

