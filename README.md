# IslandBiogeography
Assembly and annotation scripts


Assemply.IterativeBlastSearch.sh _____ This shell script batch file execute the assembly and the iterative BLASTn search. it calls Fq2Fa_interleave-pe_reads.py // fasta_length_extractor.min.max.py // fasta_id_extractor.MetaGenomicVersion.py // fasta_id_extractor.py

1.CodingFraction_dist.sh _____ This shell script uses as input file the annotation generated with MITOS WebServer and filters putative mtDNA scaffolds based on the coding fraction.

2.GeneOrientation.sh _____ This shell script uses as input scaffolds from the previous step and filters putative mtDNA scaffolds based on gene orientation.

3.GeneOrder.sh _____ This shell script uses as input scaffolds from the previous step and filters putative mtDNA scaffolds based on gene order.

4.GetAnnFasta.sh ___ This shell script uses as input scaffolds from the previous step and extract single fasta files from single genes.
