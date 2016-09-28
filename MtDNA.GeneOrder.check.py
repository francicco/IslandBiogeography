#!/usr/bin/python

from sys import argv
from collections import defaultdict


def geneOrder(geneOrder, geneOrder_dct):
	for i in range(0,len(geneOrder)):
		if i == 0:
			geneOrder_dct[geneOrder[i]].append((geneOrder[-1], geneOrder[i+1]))
		elif i == len(geneOrder)-1:
			geneOrder_dct[geneOrder[i]].append((geneOrder[i-1], geneOrder[0]))
		else:
			geneOrder_dct[geneOrder[i]].append((geneOrder[i-1], geneOrder[i+1]))
	return geneOrder_dct

# Order to check
contig=open(argv[1], 'r').readlines()

# Known order
Order=open(argv[2], 'r').readline()

#
cntg_id=argv[3]

# Strand
strand=argv[4]

Order=Order.strip().replace('-','').split()

TrueOrder_dict=defaultdict(list)

# Giving the gene order list 'argv[2]', this function create the dictionary
# that will be used to check the order of the given contig
geneOrder(Order, TrueOrder_dict)

ContigOrder=[]


#for gene in TrueOrder_dict.keys():
#	print gene, TrueOrder_dict[gene]

if len(contig) > 1:
	for gene in contig:
		gene=gene.strip()
		ContigOrder.append(gene)
#		print gene

	num_ann_genes=len(ContigOrder)

	ordered=0

	for i in range(0,len(ContigOrder)):
		gene=ContigOrder[i]
		if i == 0:
			n_gene=ContigOrder[i+1]
			if strand == '+':
				if n_gene == TrueOrder_dict[gene][0][1]:
					ordered+=1
					#print gene, n_gene, TrueOrder_dict[gene]
				#else: print gene, 'Not in the right order'
			elif strand == '-':
				if n_gene == TrueOrder_dict[gene][0][0]:
					ordered+=1
					#print gene, n_gene, TrueOrder_dict[gene]
				#else: print gene, 'Not in the right order'
		elif i == len(ContigOrder)-1:
			p_gene=ContigOrder[i-1]
			if strand == '+':
				if p_gene == TrueOrder_dict[gene][0][0]:
					ordered+=1
					#print gene, p_gene, TrueOrder_dict[gene]
				#else: print gene, 'Not in the right order'
			elif strand == '-':
				if p_gene == TrueOrder_dict[gene][0][1]:
					ordered+=1
					#print gene, p_gene, TrueOrder_dict[gene]
				#else: print gene, 'Not in the right order'
		else:
			n_gene=ContigOrder[i+1]
			p_gene=ContigOrder[i-1]
			if strand == '+':
				if (n_gene == TrueOrder_dict[gene][0][1]) and (p_gene == TrueOrder_dict[gene][0][0]):
					ordered+=1
					#print gene, p_gene, n_gene, TrueOrder_dict[gene]
				#else: print gene, 'Not in the right order'
			elif strand == '-':
				if (n_gene == TrueOrder_dict[gene][0][0]) and (p_gene == TrueOrder_dict[gene][0][1]):
					ordered+=1
					#print gene, p_gene, n_gene, TrueOrder_dict[gene]
				#else: print gene, 'Not in the right order'
	

	print '%s: %i genes, %i are ordered. Strand %s' % (cntg_id, len(ContigOrder), ordered, strand)

else:
	ContigOrder.append(contig[0].strip())
	print '%s: contains only %i gene: %s on strand %s' % (cntg_id, len(ContigOrder), ContigOrder[0], strand)
