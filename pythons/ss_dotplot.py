#!/usr/bin/env python

"""Calculate the mfes and the partition function of the sequences for both temperatures."""

# from __future__ import print_function
import sys
import RNA
import matplotlib
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

"""godzilla = Seq('GAGACCCGTAAAAGGGTCTCGAAA', IUPAC.unambiguous_dna)
foldgod = RNA.fold_compound(str(godzilla))
print dir(foldgod)"""
    
md_low = RNA.md()
md_high = RNA.md()

md_low.temperature = 30
md_high.temperature = 37

print md_low

for seq_file in SeqIO.parse(sys.stdin, "fasta"):
    # obtain the sequence string
    sequ = str(seq_file.seq)
    # generate the fold compounds for both temperatures
    fc_low = RNA.fold_compound(sequ, md_low)
    fc_high = RNA.fold_compound(sequ, md_high)
    # calculate mfe structure and mfe
    struct_low, mfe_low = fc_low.mfe()
    struct_high, mfe_high = fc_high.mfe()
    # calculate partition function for both temperatures
    pfstruct_low, pf_low = fc_low.pf()
    pfstruct_high, pf_high = fc_high.pf()
    # print sequence name, mfes and partition functions for both temperatures
    print seq_file.id
    print mfe_low, mfe_high
    print pf_low, pf_high
    

print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)
print 'matplotlib version %s' % (matplotlib.__version__)
