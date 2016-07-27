#!/usr/bin/env python

"""Generate new sequence variation: switches the first base pair of a structure"""

# from __future__ import print_function
import sys
import RNA
from Bio import SeqIO

"""Switches the first base pair of the first stem loop of a sequence"""

def determine_bp( structure ):
    # gives out the position of the first basepair in a structure
    ptable = RNA.ptable(structure)[1:]
    a, b = 0, 0
    for index, item in enumerate(ptable):
        if item:
            a = index + 1
            b = item
            break
    return a, b

def switch_bp( seqrec, a, b ):
    # switches basepair in SeqRecord indicated by the positions a and b to their complement
    complement = seqrec.reverse_complement()[::-1]
    if a < b:
        new = seqrec[:a-1] + complement[a-1] + seqrec[a:b-1] + complement.seq[b-1]+ seqrec.seq[b:]
    else:
        new = seqrec
        print 'WARNING: No valid basepairs in this structure'
    return new

for original in SeqIO.parse(sys.stdin, "fasta"):
    # obtain the sequence string
    sequ = str(original.seq)
    # generate the fold compound
    fc = RNA.fold_compound(sequ)
    # calculate mfe structure and mfe
    struct, mfe = fc.mfe()
    # change base pair
    a, b = determine_bp(struct)
    new = switch_bp(original,a,b)
    
 
    print 'old seq %s' % (original.seq)
    print 'new seq %s' % (new.seq)
    print 'struct. %s' % (struct)
    

