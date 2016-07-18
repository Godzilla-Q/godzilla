#!/usr/bin/env python

"""calculate and plot the base pair probability matrix"""

# from __future__ import print_function

import RNA
import matplotlib
import Bio
import matplotlib.pyplot as pyplot
from Bio.Seq import Seq

# Define sample sequence
godzilla = Seq('GAGACCCGTAAAAGGGTCTCGAAAGTGTGTAAAAAACACAC')
godzilla.id = 'Godzilla Queen of Monsters'
foldgod = RNA.fold_compound(str(godzilla))

# calculate partition function for both temperatures
pfstruct, pf = foldgod.pf()

# calculate base pair probability matrix
bppm = foldgod.bpp()

# print sequence name, mfes and partition functions for both temperatures
print godzilla.id

# plot base pair probability matrix
pyplot.matshow(bppm, fignum=godzilla.id)

pyplot.show()


print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)
print 'matplotlib version %s' % (matplotlib.__version__)

