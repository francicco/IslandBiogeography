#!/bin/bash


#####################################
#### Take as input a blast formated database and a fasta file with contigs/scaffolds  
#### to keep sequences matching to the database and add them to the reference DB
####
#### Usage: bash IterativeBlastSearch.sh Reads1 Reads2 MtDNA_DB Prefix Threads
####        - Reads 1 = Complete path to paired reads 1 
####        - Reads 2 = Complete path to paired reads 2 
####        - MtDNA_D = Blast formated database (should be in the same directory as the bin directory with the python scripts)  
####        - Prefix = Identifier of the library 
####        - Threads = Number of threads to use                       		
####
#### Dependencies: biopython, blastn, blastall, formatdb, fasta_id_extractor.py (bin), fasta_length_extractor.min.max.py (bin), fasta_id_extractor.MetaGenomicVersion.py (bin), IDBA-UD
####
#### Francesco Cicconardi (2016)
#####################################


### Get script arguments
R1=$1 				# Paired reads 1
R2=$2				# Paired reads 2
MTDNA_DB=$3			# mtDNA database
PREFIX=$4			# Identifier of the library/sample
THR=$5				# Number of threads

### Important parameters
MINK='23'			# Minimum kmer length for assembly
MAXK='123'			# Maximum kmer length for assembly
STEP='10'			# Increment of k-mer of each iteration
SIMILARITY='0.99'   # Similarity for alignment
MM='0'				# Max mismatch of error correction
ITERATIONS=100		# Number of iterations
E=10				# E-value for BLAST search
MIN_LEN=1000		# Minimum length of contig we want to keep [bp]
MAX_LEN=19000	    # Maximum length of contig we want to keep [bp]


### Create new directory

mkdir $PREFIX.tmp
mkdir Filtered.out
mkdir $PREFIX


### Running IDBA-UD assembly
echo
echo IDBA-UD assembly
echo Merge reads to assembly
date
fq2fa --merge $R1 $R2 $PREFIX.Trim.InterL.fasta
date
idba_ud -o $PREFIX --mink $MINK --maxk $MAXK --step $STEP -r $PREFIX.Trim.InterL.fasta --num_threads $THR --min_contig $MIN_LEN --similar $SIMILARITY --max_mismatch $MM

### Remove scaffolds/contigs according to a specified size

python bin/fasta_length_extractor.min.max.py $PREFIX/scaffold.fa $MIN_LEN $MAX_LEN > $PREFIX.tmp/$PREFIX.putative.mtDNA.fasta


### Blast search with mtDNA file as database
echo
blastall -K 1 -p blastn -d $MTDNA_DB -i $PREFIX.tmp/$PREFIX.putative.mtDNA.fasta -e $E -m 8 -o $PREFIX.tmp/$PREFIX.putative.mtDNA.blastn -n T


### Extracting the name of the contigs matching against the database
cut -f1 $PREFIX.tmp/$PREFIX.putative.mtDNA.blastn | grep -v '^#' | sort -u > $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids
python bin/fasta_id_extractor.MetaGenomicVersion.py $PREFIX.tmp/$PREFIX.putative.mtDNA.fasta $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids $MIN_LEN $MAX_LEN > $PREFIX.tmp/$PREFIX.putative.mtDNA.selected.fasta
python bin/fasta_id_extractor.py -f $PREFIX.tmp/$PREFIX.putative.mtDNA.fasta -l $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids -m exclude > $PREFIX.tmp/$PREFIX.putative.mtDNA.unselected.fasta

### Adding the matching contigs in the database
for i in $(eval echo "{1..$ITERATIONS}");
do
        if [ -f $PREFIX.tmp/$PREFIX.putative.mtDNA.selected.fasta ]; then
                mv $PREFIX.tmp/$PREFIX.putative.mtDNA.selected.fasta $PREFIX.tmp/reference_contigs.fasta
                mv $PREFIX.tmp/$PREFIX.putative.mtDNA.unselected.fasta $PREFIX.tmp/contigs.to.blast.fasta
        fi
        if [ -f $PREFIX.tmp/contigs.to.blast.selected.fasta ]; then
                mv $PREFIX.tmp/contigs.to.blast.selected.fasta $PREFIX.tmp/reference_contigs.fasta
                mv $PREFIX.tmp/contigs.to.blast.unselected.fasta $PREFIX.tmp/contigs.to.blast.fasta
        fi

        echo -e "####################################"
        echo Start Iteration ${i}
        echo

        touch $PREFIX.tmp/Iteration.${i}.started

        echo -e "\nCreating a new BLASTn Reference"
        formatdb -i $PREFIX.tmp/reference_contigs.fasta -p F -n $PREFIX.tmp/reference_contigs.mtDNA


### Blast search with new mtDNA file as database
        echo -e "\nMegaBLAST search... contigs.to.blast.fasta on reference_contigs.mtDNA"
        date
        blastall -b 1 -p blastn -d $PREFIX.tmp/reference_contigs.mtDNA -i $PREFIX.tmp/contigs.to.blast.fasta -e $E -m 8 -o $PREFIX.tmp/$PREFIX.blastn -n T

        echo -e "\nExtract aligned contigs headers"
        date
        grep -v '^#' $PREFIX.tmp/$PREFIX.blastn | cut -f1 | sort -u > $PREFIX.tmp/Iteration.${i}.contigs.headers

        cat $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids $PREFIX.tmp/Iteration.${i}.contigs.headers > $PREFIX.tmp/tmp
        mv $PREFIX.tmp/tmp $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids

        python bin/fasta_id_extractor.MetaGenomicVersion.py $PREFIX.tmp/contigs.to.blast.fasta $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids $MIN_LEN $MAX_LEN > $PREFIX.tmp/contigs.to.blast.selected.fasta
		python bin/fasta_id_extractor.py -f $PREFIX.tmp/contigs.to.blast.fasta -l $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids -m exclude > $PREFIX.tmp/contigs.to.blast.unselected.fasta


        if [ -s $PREFIX.tmp/Iteration.${i}.contigs.headers ]; then
                echo -e "Start Iteration ${i} ended"
                touch $PREFIX.tmp/Iteration.${i}.ended
        else
                echo -e "\nNo more contigs to align"
                touch $PREFIX.tmp/Iterations.ended
				python bin/fasta_id_extractor.py -f $PREFIX.tmp/$PREFIX.putative.mtDNA.fasta -l $PREFIX.tmp/$PREFIX.ColNTs.blastn.ids -m select > Filtered.out/$PREFIX.filtered.mtDNA.contigs.fa
				break
        fi


done







