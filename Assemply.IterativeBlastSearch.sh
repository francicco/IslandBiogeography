#!/bin/bash


WORKING_DIR=$1
TMP_DIR=/home/lv70640/c7701100/.tmp_dir
TAG='IDBA-UD.All.assembly'
READS1='merged.Azor_Cana.30smpls_R1_001'
READS2='merged.Azor_Cana.30smpls_R2_001'
MTDNA_DB=$1/Annotation/Arthropoda.mtDNA
THREADS=$2
MEM=$3

MINK='23'
MAXK='123'
STEP='10'
MINCONTIG='1000'
MAXCONTIG='19000'
SIMILARITY='0.99'
MM='0'
E=10 #e value for BLAST search
ITERATIONS=100

cd $WORKING_DIR

mkdir -p $TAG

echo
echo IDBA-UD assembly
echo Merge reads to assembly
date
Fq2Fa_interleave-pe_reads.py $READS1.fastq $READS2.fastq > $TAG/pe-reads.fasta
date
idba_ud -o $TMP_DIR/$TAG --mink $MINK --maxk $MAXK --step $STEP -r $TAG/pe-reads.fasta --num_threads $THREADS --min_contig $MINCONTIG --similar $SIMILARITY --max_mismatch $MM

cp -r $TMP_DIR/$TAG .
rm -r $TMP_DIR/$TAG


cd $TAG


echo
echo MegaBLAST search... reads.fasta on $MTDNA_DB2
date
fasta_length_extractor.min.max.py scaffold.fa $MINCONTIG $MAXCONTIG > $TAG.putative.mtDNA.fasta
blastall -K 1 -P 0 -p blastn -d $MTDNA_DB -i $TAG.putative.mtDNA.fasta -e $E -m 8 -o $TAG.putative.mtDNA.blastn -n T

cut -f1 $TAG.ColNTs.blastn | grep -v '^#' | sort -u > $TAG.ColNTs.blastn.ids
fasta_id_extractor.MetaGenomicVersion.py $TAG.putative.mtDNA.fasta $TAG.ColNTs.blastn.ids

for i in $(seq $ITERATIONS)
do
        if [ -f $TAG.putative.mtDNA.selected.fasta ]; then
                mv $TAG.putative.mtDNA.selected.fasta reference_contigs.fasta
                mv $TAG.putative.mtDNA.unselected.fasta contigs.to.blast.fasta
        fi
        if [ -f contigs.to.blast.selected.fasta ]; then
                mv contigs.to.blast.selected.fasta reference_contigs.fasta
                mv contigs.to.blast.unselected.fasta contigs.to.blast.fasta
        fi

        echo -e "####################################"
        echo Start Iteration ${i}
        echo

        touch Iteration.${i}.started

        echo -e "\nCreating a new BLASTn Reference"
        formatdb -i reference_contigs.fasta -p F -n reference_contigs.mtDNA

        echo -e "\nMegaBLAST search... contigs.to.blast.fasta on reference_contigs.mtDNA"
        date
        blastall -a $THREADS -b 1 -p blastn -d reference_contigs.mtDNA -i contigs.to.blast.fasta -e $E -m 8 -o $TAG.blastn -n T

        echo -e "\nExtract alined contigs headers"
        date
        grep -v '^#' $TAG.blastn | cut -f1 | sort -u > Iteration.${i}.contigs.headers

        cat $TAG.ColNTs.blastn.ids Iteration.${i}.contigs.headers > tmp
        mv tmp $TAG.ColNTs.blastn.ids

        fasta_id_extractor.MetaGenomicVersion.py contigs.to.blast.fasta $TAG.ColNTs.blastn.ids

        if [ -s Iteration.${i}.contigs.headers ]; then
                echo -e "Start Iteration ${i} ended"
                touch Iteration.${i}.ended
        else
                echo -e "\nNo more contigs to align"
                touch Iterations.ended
                break
        fi
done


fasta_id_extractor.py $TAG.putative.mtDNA.fasta $TAG.ColNTs.blastn.ids > IDBA-UD.qft.mtDNA.scaffolds.fa










