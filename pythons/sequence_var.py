#!/usr/bin/env python

"""
Convert secondary structure to pair table.
generate new sequence variation
    convert pair table (array) to my own secondary structure file
    switch nts of one basepair
    
    """

# from __future__ import print_function
import sys
import RNA
import matplotlib
from math import exp
import re # for regular expressions
import matplotlib.pyplot as plt
from math import exp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from matplotlib.backends.backend_pdf import PdfPages

"""Switches the first base pair of the first stem loop of a sequence"""


def switch_bp( seqrec, a, b ):
    complement = seqrec.reverse_complement()[::-1]
    new = seqrec[:a-1] + complement[a-1] + seqrec[a:b-1] + complement.seq[b-1]+ seqrec.seq[b:]
    return new    

for original in SeqIO.parse(sys.stdin, "fasta"):
    # obtain the sequence string
    sequ = str(original.seq)
    # generate the fold compounds
    fc = RNA.fold_compound(sequ)
    # calculate mfe structure and mfe
    struct, mfe = fc.mfe()
    # calculate partition function
    pfstruct, pf = fc.pf()
    # calculate base pair probability matrix
    bppm = fc.bpp()
    # change base pair
    ptable = RNA.ptable(struct)[1:]
    a, b = 0, 0
    for index, item in enumerate(ptable):
        if item:
            a = index + 1
            b = item
            break
    new = switch_bp(original,a,b)
 
    print original.seq
    print new.seq
    print struct
    

