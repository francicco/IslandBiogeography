#!/usr/bin/python

from Bio import SeqIO
from sys import argv
from collections import defaultdict



fasta_file=open(argv[1], 'r')

min_length=int(argv[2])
max_length=int(argv[3])

for record in SeqIO.parse(fasta_file, 'fasta'):
	if len(record.seq) >= min_length and len(record.seq) <= max_length:
		print '>'+record.description
		print record.seq
