#!/usr/bin/python

import optparse
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import defaultdict


################################# Command line options

desc='Generate the consensus of an alignment.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 31-03-2015 - Author: FCicconardi')

parser.add_option('-i', '--alignment', dest='aln', help='Input alignment in FASTA. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-t', '--consensus-threshold', dest='thld', help='Set consensus threshold value. The higher the more is sensitive [0-1]. Default <0.6>.', action='store', metavar='<ARG>', type='float', default=0.4)

(opts, args) = parser.parse_args()

mandatories = ['aln']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING!  output file is not specified\n"
                parser.print_help()
                exit(-1)


############################## Reading files and parameters

alignment = AlignIO.read(open(opts.aln), 'fasta')
print "Alignment length %i" % alignment.get_alignment_length()

summary_align = AlignInfo.SummaryInfo(alignment)

consensus = str(summary_align.dumb_consensus(threshold=float(opts.thld), ambiguous='-', require_multiple = 2)).replace('!','-')

consensus=''
for i in range(0,alignment.get_alignment_length()):
	pos=''
	amb=[]
	for record in alignment:
		pos+=record.seq[i].replace('!','-')
	pA=float(pos.count('A'))/len(pos)
	pT=float(pos.count('T'))/len(pos)
	pC=float(pos.count('C'))/len(pos)
	pG=float(pos.count('G'))/len(pos)
	pGap=float(pos.count('-'))/len(pos)
	#print pos
	#print i, 'A', pA
	#print i, 'T', pT
	#print i, 'C', pC
	#print i, 'G', pG
	#print i, '-', pGap
	if pA > opts.thld:
		#print i, 'A', pA
		amb.append('A')
	if pT > opts.thld:
		#print i, 'T', pT
		amb.append('T') 
	if pC > opts.thld:
		#print i, 'C', pC
		amb.append('C')
	if pG > opts.thld:
		#print i, 'G', pG
		amb.append('G')
	if pGap > opts.thld:
		#print i, '-', pGap
		amb.append('-')
	if len(amb) == 0 or len(amb) == 5:
		consensus+='N'
	elif len(amb) == 1:
		consensus+=amb[0]
	elif len(amb) == 2:
		if 'A' in amb and 'G' in amb:
			consensus+='R'
		elif 'C' in amb and 'T' in amb:
			consensus+='Y'
		elif 'G' in amb and 'C' in amb:
			consensus+='S'
		elif 'A' in amb and 'T' in amb:
			consensus+='W'
		elif 'T' in amb and 'G' in amb:
			consensus+='K'
		elif 'A' in amb and 'C' in amb:
			consensus+='M'
		elif 'A' in amb and '-' in amb:
			consensus+='A'
		elif 'T' in amb and '-' in amb:
			consensus+='T'
		elif 'C' in amb and '-' in amb:
			consensus+='C'
		elif 'G' in amb and '-' in amb:
			consensus+='G'
	elif len(amb) == 3:
		if 'C' in amb and 'G' in amb and 'T' in amb:
			consensus+='B'
		elif 'A' in amb and 'G' in amb and 'T' in amb:
			consensus+='D'
		elif 'A' in amb and 'C' in amb and 'T' in amb:
			consensus+='H'
		elif 'A' in amb and 'C' in amb and 'G' in amb:
			consensus+='V'
		elif 'A' in amb and 'T' in amb and '-' in amb:
			consensus+='W'
		elif 'A' in amb and 'C' in amb and '-' in amb:
			consensus+='M'
		elif 'A' in amb and 'G' in amb and '-' in amb:
			consensus+='R'
		elif 'T' in amb and 'C' in amb and '-' in amb:
			consensus+='Y'
		elif 'T' in amb and 'G' in amb and '-' in amb:
			consensus+='K'
		elif 'C' in amb and 'G' in amb and '-' in amb:
			consensus+='S'
	elif len(amb) == 4:
		if '-' not in amb:
			consensus+='N'
		elif 'A' in amb and 'T' in amb and 'C' in amb and '-' in amb:
			consensus+='H'
		elif 'A' in amb and 'C' in amb and 'G' in amb and '-' in amb:
			consensus+='V'
		elif 'A' in amb and 'T' in amb and 'G' in amb and '-' in amb:
			consensus+='D'
		elif 'T' in amb and 'C' in amb and 'G' in amb and '-' in amb:
			consensus+='B'


tmp=opts.aln.split('.alned')
locus=tmp[0]

filename=locus+'.consensus.fasta'

out_fasta=open(filename, 'w')

print >> out_fasta, '>'+opts.aln.replace('.fasta', '')+'_consensus threshold '+str(opts.thld)
print >> out_fasta, consensus



#print opts.thld


	
