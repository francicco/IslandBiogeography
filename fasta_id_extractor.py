#!/usr/bin/python

import optparse
from sys import argv
from Bio import SeqIO
from collections import defaultdict

################################# Command line options

desc='Fasta filterer'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 09-03-2015 - Author: FCicconardi')

parser.add_option('-f', '--fasta-file', dest='fasta', help='Fasta file you want to filter. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-l', '--id-list', dest='ids', help='Id list to filter in or out. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-m', '--method', dest='mtd', help='Select or Exclude those ids from the fasta file (select || exclude). Mandatory opt.', action='store', metavar='<ARG>')


(opts, args) = parser.parse_args()

mandatories = ['fasta','ids', 'mtd']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING! One or more options not specified\n"
                parser.print_help()
                exit(-1)

############################## Reading files and parameters

fasta_file=open(opts.fasta, 'r')

id_list=open(opts.ids, 'r').readlines()
idlist=[]

for el in id_list:
	idlist.append(el.strip())

for record in SeqIO.parse(fasta_file, 'fasta'):
	if opts.mtd == 'select':
		if record.id in idlist:
			print '>'+record.description
			print record.seq
	elif opts.mtd == 'exclude':
		if record.id in idlist: continue
		else:
			print '>'+record.description
			print record.seq



