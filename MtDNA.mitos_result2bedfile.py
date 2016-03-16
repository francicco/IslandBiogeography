#!/usr/bin/python

from sys import argv
from collections import defaultdict

result=open(argv[1], 'r').readlines()


# conversion
for row in result:
	row=row.strip().split('\t')
	scf=row[0]
	score=row[7]
	if row[6] == '1': strand='+'
	elif row[6] == '-1': strand='-'
	start=row[4]
	end=row[5]
	name=row[2]
	print '%s\t%s\t%s\t%s\t%s\t%s' % (scf,start,end,name,score,strand)
		



