#!/usr/bin/env python

"""plot the partition function"""

# from __future__ import print_function
import sys
import RNA
import matplotlib
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

godzilla = Seq('GAGACCCGTAAAAGGGTCTCGAAA', IUPAC.unambiguous_dna)
godzilla.id = 'Godzilla Queen of Monsters'
monster = str(godzilla)
foldgod = RNA.fold_compound(monster)

# calculate mfe structure and mfe
struct, mfe = foldgod.mfe()

# calculate partition function for both temperatures
pfstruct, pf = foldgod.pf()

    # print sequence name, mfes and partition functions for both temperatures
print godzilla.id
print mfe
print pf
    

print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)
print 'matplotlib version %s' % (matplotlib.__version__)
