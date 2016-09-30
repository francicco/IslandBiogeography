#!/usr/bin/python

import optparse
from os import listdir
from subprocess import call
from Bio import SeqIO, AlignIO
from collections import defaultdict

################################# Command line options

desc='Convert alignment from FASTA to NEXUS format.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 31-03-2015 - Author: FCicconardi')
parser.add_option('-i', '--input-directory', dest='dir', help='Directory where all alignments are (FASTA format). Default current directory.', action='store', metavar='<FILE>', default='.')
parser.add_option('-s', '--suffix', dest='sfx', help='Suffix of alignments', action='store', metavar='<ARG>')
parser.add_option('-o', '--nexus-output', dest='nex', help='Output NEXUS file.', action='store', metavar='<FILE>')

(opts, args) = parser.parse_args()

mandatories = ['sfx','nex']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING!  output file is not specified\n"
                parser.print_help()
                exit(-1)


############################## Reading files and parameters


## Create a list with all fasta files
directory=listdir(opts.dir)

files=[]

for file in directory:
	if file.endswith(opts.sfx):
		files.append(file)

#####################################


## Create a set with all contigs name contained in all fasta file

samples_dict=defaultdict(list)

taxa=[]

for file in files:
	print file
	for record in SeqIO.parse(open(file, 'r'), 'fasta'):
		id=record.id.replace(';','')
		taxa.append(id)

taxaList=set(taxa)

ntaxa=len(taxaList)

#################################################################

## FASTA conversion

all_seq_loci=defaultdict(list)

for file in files:
	tmp_dict=defaultdict(list)
	for record in SeqIO.parse(open(file, 'r'), 'fasta'):
		id=record.id
		tmp=record.description.split()
		tmp_dict[id].append(str(record.seq).replace('!','-'))
	tmp_=file.split('.')
	#locus=tmp[-1]
	locus=tmp_[2]
	loc_tax=len(tmp_dict.keys())
	nchar=len(record.seq)
	mis_len='?'*nchar
	tmp_out=open('.'+locus+'.'+str(nchar)+'.seqs', 'w')
	output_nexus=open(file.replace('.fasta','.nex'), 'w')
	print 'Processing locus %s, %i taxa' % (locus, loc_tax)
	print >> output_nexus, '''#NEXUS
	begin data;
	dimensions ntax=%s nchar=%i;
	format datatype=dna missing=? gap=-;
	matrix''' % (loc_tax, nchar)
	for id in taxaList:
                if id in tmp_dict.keys():
                        print >> tmp_out, id+'\t'+tmp_dict[id][0]
			print >> output_nexus, id+'\t'+tmp_dict[id][0]
		else:
                        print >> tmp_out, id+'\t'+mis_len
	print >> output_nexus, ';\nend;'
	print >> output_nexus, '''begin sets;
charset %s-1 = 1-%i\\3;
charset %s-2 = 2-%i\\3;
charset %s-3 = 3-%i\\3;
charpartition combined = p1: %s-1, p2: %s-2, p3: %s-3;
end;''' % (locus, nchar, locus, nchar, locus, nchar, locus, locus, locus)
	tmp_out.close()
	output_nexus.close()		

###############################################################################

## Concatenation

## Create a dictionary with all tmp files to concatenate

directory=listdir(opts.dir)

all_dict=defaultdict(list)

for file in directory:
        if file.endswith('.seqs'):
		tmp=file.split('.')
		locus=tmp[1]
		nchar=tmp[2]
		file=open(file, 'r').readlines()
		for row in file:
			row=row.strip().split('\t')
			all_dict[row[0]].append((locus+' '+nchar, row[1]))
		call(['rm', '.'+locus+'.'+nchar+'.seqs'])

#####################################

concat=[]

ntaxa=len(all_dict.keys())

for scf in all_dict.keys():
	locus_order=''
	seqs_conc=''
	for l in range(0,len(all_dict[scf])):
		locus_order+=all_dict[scf][l][0]+'; '
		seqs_conc+=all_dict[scf][l][1]
	nchar=len(seqs_conc)
	concat.append(scf+'\t'+seqs_conc)
	

allseqs=''


conc_nexus_out=open(opts.nex, 'w')

print >> conc_nexus_out, '''#NEXUS
begin data;
dimensions ntax=%s nchar=%i;
format datatype=dna missing=? gap=-;
matrix''' % (ntaxa, nchar)
for id in concat:
	seq=id.split('\t')
	allseqs+=seq[1]
	print >> conc_nexus_out, id
print >> conc_nexus_out, ';\nend;\n\nbegin sets;'

loci=locus_order[:-2].split(';')

end=0
char_comb=''
c=0
for el in loci:
	el=el.split()
	locus=el[0]
	lchar=int(el[1])
	if locus.startswith('rrn'):
		start=end+1
		end=end+lchar
		print >> conc_nexus_out, 'charset %s = %i-%i;' % (locus, start, end)
		c+=1
		set='p%i: %s, ' % (c, locus)
		char_comb+=set
	else:
		fpos=end+1
		spos=end+2
		tpos=end+3
		end=end+lchar
		print >> conc_nexus_out, 'charset %s-1  = %i-%i\\3;' % (locus, fpos, end)
		c+=1
		set1='p%i: %s-1, ' % (c, locus)
		print >> conc_nexus_out, 'charset %s-2  = %i-%i\\3;' % (locus, spos, end)
		c+=1
		set2='p%i: %s-2, ' % (c, locus)
		print >> conc_nexus_out, 'charset %s-3  = %i-%i\\3;' % (locus, tpos, end)
		c+=1
		set3='p%i: %s-3, ' % (c, locus)
		char_comb+=set1+set2+set3
print >> conc_nexus_out, 'charpartition combined = '+char_comb[:-2]+';\nend;'

conc_nexus_out.close()

gaps=allseqs.count('-')
missing=allseqs.count('?')
allchar=len(allseqs)
nucl=allchar-gaps-missing

print '''\n\n\n%i, Taxa for %i characters:\n
Total amount of characters %i (100%%), of them:
%i (%g%%) are defined nucleotides
%i (%g%%) are gaps, and
%i (%g%%) are missing.''' % (ntaxa, nchar, allchar, nucl, float(nucl)/float(allchar)*100, gaps, float(gaps)/float(allchar)*100, missing, float(missing)/float(allchar)*100)

