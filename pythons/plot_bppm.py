#!/usr/bin/env python

"""calculate and plot the base pair probability matrix"""

# from __future__ import print_function

import RNA
import matplotlib
import Bio
import matplotlib.pyplot as plt
from Bio.Seq import Seq

# Define sample sequence
godzilla = Seq('GAGACCCGTAAAAGGGTCTCGAAAGTGTGTAAAAAACACAC')
godzilla.id = 'Godzilla Queen of Monsters'
#foldgod = RNA.fold_compound(str(godzilla))
Hirsch = Seq('CCGCACAGCGGGCAGUGCCC')
Hirsch.id = 'Papa Hirsch protects us all'
#foldHirsch = RNA.fold_compound(str(Hirsch))

monsters = (godzilla, Hirsch)
# use either 'BuPu' or 'Greys'
colormap = 'BuPu'

def plot_bppm ( bppm, name ):
# plot base pair probability matrix, write plot to post script file
    plt.matshow(bppm, fignum=name, cmap=plt.get_cmap(colormap))
    plt.savefig('%s.ps' % (name), format='ps')
    plt.close()
    return

for monster in monsters:
# calculate 1) foldcompound 2) partition function 3) base pair probability matrix in that order (!)
    foldmonster = RNA.fold_compound(str(monster))
    pfstruct, pf = foldmonster.pf()
    bppm = foldmonster.bpp()
    plot_bppm(bppm, monster.id)

print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)
print 'matplotlib version %s' % (matplotlib.__version__)

