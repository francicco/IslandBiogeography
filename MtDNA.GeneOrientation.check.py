#!/usr/bin/python

from sys import argv
from collections import defaultdict

orientation=open(argv[1], 'r').readline().strip().split(' ')

bad_file=open(argv[2], 'r').readlines()

gene_dict=defaultdict(list)

fstCycl=defaultdict(list)

cntg_g_order=open('cntg_g_order.dat', 'w')

tmp_genes_dict=defaultdict(list)

for el in bad_file:
    rcd=el.strip().split('\t')
    qval=float(rcd[7])
    if rcd[2].startswith('trn'): continue
    else:
        if rcd[6] == '1':
            print >> cntg_g_order, rcd[2]
            fstCycl[rcd[2]].append((rcd[7], int(rcd[4]), int(rcd[5]), '+'))
        else:
            print >> cntg_g_order, rcd[2]
            fstCycl[rcd[2]].append((rcd[7], int(rcd[5]), int(rcd[4]), '-'))




cntg_g_order.close()

c=0

for gene in orientation:
    locus=gene[:-1]
    strand=gene[-1]
    if locus in fstCycl.keys():
        if strand == fstCycl[locus][0][3]:
            c+=1


if c == 0:
    #print '100%% are consistenty oriented, but the opposite strand. %i out of %i.' % (c, len(sndCycl.keys()))
    print argv[3]+' -'
elif c == len(fstCycl.keys()):
#    #print '100%% are consistenty oriented, %d out of %d.' % (c, len(sndCycl.keys()))
    print argv[3]+' +'
else:
    match=float(c)
    total=float(len(fstCycl.keys()))
#    #print 'Not consistet orientation of %d out of %d.' % (c, len(sndCycl.keys()))
