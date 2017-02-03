#!/usr/bin/python

from Bio import SeqIO
from sys import argv
from collections import defaultdict



fasta_file=open(argv[1], 'r')

id_list=open(argv[2], 'r').readlines()
idlist=[]


for el in id_list:
	idlist.append(el.strip())

for record in SeqIO.parse(fasta_file, 'fasta'):
	if record.id in idlist:
				print '>'+record.description
				print record.seq