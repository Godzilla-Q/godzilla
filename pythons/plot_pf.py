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

# calculate base pair probability matrix
bppm = foldgod.bpp()

# print sequence name, mfes and partition functions for both temperatures
print godzilla.id
print mfe
print pf
print bppm

# plot base pair probability matrix

import matplotlib.pyplot as pyplot
import numpy

def samplemat(dims):
    """Make a matrix with all zeros and increasing elements on the diagonal"""
    aa = numpy.zeros(dims)
    for i in range(min(dims)):
        aa[i, i] = i
    return aa

# Display 2 matrices of different sizes
dimlist = [(15, 20)]
for d in dimlist:
    pyplot.matshow(samplemat(d))

# Display a random matrix with a specified figure number and a grayscale
# colormap
pyplot.matshow(numpy.random.rand(64, 64), fignum=100, cmap=pyplot.cm.gray)

pyplot.show()


print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)
print 'matplotlib version %s' % (matplotlib.__version__)

