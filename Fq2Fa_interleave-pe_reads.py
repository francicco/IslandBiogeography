#!/usr/bin/python

from sys import argv, stdout

reads1=open(argv[1], 'r')

idline1=reads1.readline()
seq1   =reads1.readline()
spacer1=reads1.readline()
quals1 =reads1.readline()

reads2=open(argv[2], 'r')

idline2=reads2.readline()
seq2   =reads2.readline()
spacer2=reads2.readline()
quals2 =reads2.readline()

while idline1 and idline2:
	idline1=idline1.replace('@','>')
	idline2=idline2.replace('@','>')
        print '%s%s' % ( idline1, seq1.strip())
	print '%s%s' % ( idline2, seq2.strip())
	idline1=reads1.readline()
        seq1   =reads1.readline()
        spacer1=reads1.readline()
        quals1 =reads1.readline()
        idline2=reads2.readline()
        seq2   =reads2.readline()
        spacer2=reads2.readline()
        quals2 =reads2.readline()


