#!/usr/bin/python

import optparse
from Bio import SeqIO
from collections import defaultdict


################################# Command line options

desc='Extract DNA sequences into a fasta file based on feature circular coordinates.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 17-03-2015 - Author: FCicconardi')

parser.add_option('-f', '--fasta-file', dest='fasta', help='Input FASTA file. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-b', '--bed-file', dest='bed', help='BED file of ranges to extract from -f. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-o', '--output-file', dest='out', help='Output fasta file. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['fasta', 'bed', 'out']
for m in mandatories:
	if not opts.__dict__[m]:
		print "\nWARNING!  output file is not specified\n"
		parser.print_help()
		exit(-1)


############################## Reading files and parameters

bed=open(opts.bed, 'r').readlines()

bed_dct=defaultdict(list)

for ann in bed:
	ann=ann.strip().split('\t')
	bed_dct[ann[3]].append((ann[1],ann[2],ann[5]))


fasta=open(opts.fasta, 'r')

out=open(opts.out, 'w')

for record in SeqIO.parse(fasta, 'fasta'):
	for ann in bed_dct:
		start=int(bed_dct[ann][0][0])
		end=int(bed_dct[ann][0][1])
		strand=bed_dct[ann][0][2]
		name='%s; %i-%i; %s; %s' % (record.id, start, end, strand, ann)
		if start < 0:
			if len(record.seq) > 12000:
				seq=record.seq[start:]+record.seq[:end]
			else:
				seq=record.seq[:end]
				name='%s; %i-%i; %s; %s' % (record.id, 0, end, strand, ann)
			if strand == '+':
				print >> out, '>'+name+'\n'+seq
			elif strand == '-':
				print >> out,  '>'+name+'\n'+seq.reverse_complement()
		else:
			if strand == '+':
				seq=record.seq[start:end]
				print >> out, '>'+name+'\n'+seq
			elif strand == '-':
				seq=record.seq[start+3:end-2]
				print >> out, '>'+name+'\n'+seq.reverse_complement()


out.close()



