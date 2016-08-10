#!/usr/bin/env python

"""calculate and plot the base pair probability matrix"""

# from __future__ import print_function

import RNA
import matplotlib
import Bio
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from ss_dotplot import versions_used

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
    plt.savefig('{:s}.ps'.format(name), format='ps')
    plt.close()
    return

def main():
    for monster in monsters:
        # calculate 1) foldcompound 2) partition function 3) base pair probability matrix in that order (!)
        foldmonster = RNA.fold_compound(str(monster))
        pfstruct, pf = foldmonster.pf()
        bppm = foldmonster.bpp()
        plot_bppm(bppm, monster.id)
    print versions_used()

if __name__ == '__main__':
    main()
