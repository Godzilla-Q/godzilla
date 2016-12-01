#!/usr/bin/env python

"""Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""

import sys
import RNA
import argparse
from Bio.SeqUtils import GC
from subprocess import Popen, PIPE
from seq_properties import argument_parser
from seq_properties import temperature_reactivity

def bp_stats( sequence, structure ): # work in progress!
    """Calculate GC content and other stats"""
    # read sequence and structure and turn them into lists I can use
    ptable = list(RNA.ptable(structure))
    length = ptable.pop(0)
    seq_l = list(sequence)
    # make anoher list containing the pairing partners
    pair_partners = []
    for i in ptable:
        if i != 0:
            pair_partners.append(sequence[i - 1])
        else:
            pair_partners.append('-')
    # make a list of tuples containing the basepairs
    pairlist = zip(seq_l, pair_partners)
    # count all possible basepairs and all unpaired bases
    AUpairs = pairlist.count(('A', 'U')) + pairlist.count(('U', 'A'))
    GUpairs = pairlist.count(('G', 'U')) + pairlist.count(('U', 'G'))
    GCpairs = pairlist.count(('G', 'C')) + pairlist.count(('C', 'G'))    
    Aunp = pairlist.count(('A', '-'))
    Cunp = pairlist.count(('C', '-'))
    Gunp = pairlist.count(('G', '-'))
    Uunp = pairlist.count(('U', '-'))
    results = [AUpairs, GUpairs, GCpairs, Aunp, Cunp, Gunp, Uunp]
    results = [float(i) / length for i in results]
    results.append(GC(sequence))
    return results

def main(): 
    """Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""
    print 'sequence, structure, energy_of_struct1, energy_of_struct2, gc_content, AU, GU, GC, Aunp, Cunp, Gunp, Uunp, delta_energy'
    parser = argument_parser(sys.argv[1:])
    temperature_1, temperature_2 = parser.temperatures
    for line in sys.stdin:
        sequence, structure = line.split()
        delta_energy, energy_of_struct1, energy_of_struct2 = temperature_reactivity(sequence, structure, temperature_1, temperature_2)
        [AU, GU, GC, Aunp, Cunp, Gunp, Uunp, GCcontent] = bp_stats(sequence, structure)
        print sequence, structure, energy_of_struct1, energy_of_struct2, GCcontent, AU, GU, GC, Aunp, Cunp, Gunp, Uunp, delta_energy

    
if __name__ == '__main__':
    sys.exit(main())
