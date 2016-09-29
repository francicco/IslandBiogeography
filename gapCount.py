#!/usr/bin/python

from sys import argv
from Bio import SeqIO


fasta=open(argv[1], 'r')
gap_pec=float(argv[2])

#fst_out=open(argv[3], 'w')

for record in SeqIO.parse(fasta, 'fasta'):
	gaps=float(record.seq.count('-'))/float(len(record.seq))
	if gaps <= gap_pec:
		#print >> fst_out, '>%s\n%s' % (record.description, str(record.seq).replace('-',''))
		print '>%s\n%s' % (record.description, str(record.seq).replace('-',''))
	else:
		print record.description, gaps

#fst_out.close() 
