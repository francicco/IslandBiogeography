#!/usr/bin/python

from Bio import SeqIO
from sys import argv
from collections import defaultdict



fasta_file=open(argv[1], 'r')

id_list=open(argv[2], 'r').readlines()
idlist=[]

out_fasta=open(argv[3], 'w')

for el in id_list:
	idlist.append(el.strip())

for record in SeqIO.parse(fasta_file, 'fasta'):
	if record.id in idlist:
		print >> out_fasta, '>'+record.description
		print >> out_fasta, record.seq
	else:
		print '>'+record.description
