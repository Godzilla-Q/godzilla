#!/usr/bin/env python

"""Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""

import sys
import RNA
import argparse
from Bio.SeqUtils import GC
from subprocess import Popen, PIPE
from seq_properties import argument_parser
from seq_properties import temperature_reactivity

def bp_stats_old( sequence, structure ):
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
    
    # make a list indicating whether the bp/unpaired base sits in the middle or on the edge
    mrlist = []
    
    structure_l = list(structure)
    if structure_l[0] == structure_l[1]:
        mrlist.append('m')
    else:
        mrlist.append('r')
    for i, j, k in zip(structure_l, structure_l[1:], structure_l[2:]):
        if i == j == k:
            mrlist.append('m')
        else:
            mrlist.append('r')
    if structure_l[-1] == structure_l[-2]:
        mrlist.append('m')
    else:
        mrlist.append('r')
    
    mrlist.append('m')
    # split pairlist in two lists, one containing the middle pairs, one the edge pairs
    list_m = []
    list_r = []
    for i, j in zip(pairlist, mrlist):
        if j == 'm':
            list_m.append(i)
        if j == 'r':
            list_r.append(i)
    
    # count all possible basepairs and all unpaired bases in the middle list
    AU_m = list_m.count(('A', 'U')) + list_m.count(('U', 'A'))
    GU_m = list_m.count(('G', 'U')) + list_m.count(('U', 'G'))
    GC_m = list_m.count(('G', 'C')) + list_m.count(('C', 'G'))    
    A_m = list_m.count(('A', '-'))
    C_m = list_m.count(('C', '-'))
    G_m = list_m.count(('G', '-'))
    U_m = list_m.count(('U', '-'))
    results_m = [AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m]
    # count all possible basepairs and all unpaired bases in the edge list    
    AU_r = list_r.count(('A', 'U')) + list_r.count(('U', 'A'))
    GU_r = list_r.count(('G', 'U')) + list_r.count(('U', 'G'))
    GC_r = list_r.count(('G', 'C')) + list_r.count(('C', 'G'))    
    A_r = list_r.count(('A', '-'))
    C_r = list_r.count(('C', '-'))
    G_r = list_r.count(('G', '-'))
    U_r = list_r.count(('U', '-'))
    results_r = [AU_r, GU_r, GC_r, A_r, C_r, G_r, U_r]
    
    results = results_m + results_r
    results = [float(i) / length for i in results]
    results.append(GC(sequence))
    return results

def main(): 
    """Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""
    print 'sequence, structure, energy_of_struct1, energy_of_struct2, AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m, AU_r, GU_r, GC_r, A_r, C_r, G_r, U_r, GCcontent, delta_energy'
    parser = argument_parser(sys.argv[1:])
    temperature_1, temperature_2 = parser.temperatures
    for line in sys.stdin:
        sequence, structure = line.split()
        delta_energy, energy_of_struct1, energy_of_struct2 = temperature_reactivity(sequence, structure, temperature_1, temperature_2)
        [AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m, AU_r, GU_r, GC_r, A_r, C_r, G_r, U_r, GCcontent] = bp_stats(sequence, structure)
        results = [sequence, structure, energy_of_struct1, energy_of_struct2, AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m, AU_r, GU_r, GC_r, A_r, C_r, G_r, U_r, GCcontent, delta_energy]
        print ', '.join(map(str, results))
            
    
if __name__ == '__main__':
    sys.exit(main())

