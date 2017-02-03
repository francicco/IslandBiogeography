#!/usr/bin/python

import optparse
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from collections import defaultdict

def ungap(seq, gaps):
    c=0
    ungapped=''
    for n in seq:
            c+=1
            if c not in gaps:
                    ungapped+=n.replace('!','-')
    return ungapped


################################# Command line options

desc='Generate the consensus of an alignment.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 31-03-2015 - Author: FCicconardi')

parser.add_option('-i', '--alignment', dest='aln', help='Input alignment in FASTA. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-c', '--consensus', dest='cons', help='Consensus in FASTA format.', action='store', metavar='FILE')
parser.add_option('-t', '--type', dest='typ', help='Type of DNA. It could be <coding> or <noncoding>. Default: coding.', action='store', metavar='<ARG>', default='coding')
parser.add_option('-o', '--output', dest='out', help='Output ungapped alignment in FASTA format.', action='store', metavar='FILE') 

(opts, args) = parser.parse_args()

mandatories = ['aln', 'cons', 'out']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING!  output file is not specified\n"
                parser.print_help()
                exit(-1)


############################## Reading files and parameters

for record in SeqIO.parse(open(opts.cons), 'fasta'):
    c=0
    cod=1
    codon_dict=defaultdict(list)
    consensus=record.seq
    for pos in consensus:
        c+=1
        if c%3 == 1: cod+=1
        if pos == '-':
            codon_dict[cod].append((c))

gaps=[]
for cod in codon_dict.keys():
    if opts.typ=='coding':
        if len(codon_dict[cod]) > 0:
            gaps.append((cod-1)*3-2)
            gaps.append((cod-1)*3-1)
            gaps.append((cod-1)*3)
    elif opts.typ=='noncoding':
        for pos in range(0,len(codon_dict[cod])):
            gaps.append(codon_dict[cod][pos])
    
alignment = AlignIO.read(open(opts.aln), 'fasta')
print '\nAlignment length %i' % alignment.get_alignment_length()

ungapped_fasta=open(opts.out, 'w')

for record in alignment:
    ungapped=ungap(record.seq, gaps)
    print >> ungapped_fasta, '>'+record.description.replace(';','')
    print >> ungapped_fasta, ungapped 

ungapped_fasta.close()

print '\n%i Position removed.\n' % (len(gaps))
