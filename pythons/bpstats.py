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
    
    # split pairlist in lists, depending on whether the bp/unpaired base sits in the middle or on the edge and whether the unpaired base is 3' or 5' of the next basepair

    bp_m = []
    up_m = []
    bp_r = []
    up_r5 = []
    up_r3 = []
    structure_l = list(structure)
    if structure_l[0] == structure_l[1]:
        if structure[0] == '.':
            up_m.append(pairlist[0])
        else:
            bp_r.append(pairlist[0])
    else:
        if structure[0] == '.':
            up_r5.append(pairlist[0])
        else:
            bp_r.append(pairlist[0])

    for i, j, k, l in zip(structure_l, structure_l[1:], structure_l[2:], pairlist[1:]):
        if i == j == k:
            if j == '.':
                up_m.append(l)
            else:
                bp_m.append(l)
        else:
            if j == '.':
                if i == ')' or i == '(':
                    up_r3.append(l)
                if k == ')' or k == '(':
                    up_r5.append(l)
            else:
                bp_r.append(l)

    if structure_l[-1] == structure_l[-2]:
        if structure[0] == '.':
            up_m.append(pairlist[-1])
        else:
            bp_r.append(pairlist[-1])
    else:
        if structure[0] == '.':
            up_r3.append(pairlist[-1])
        else:
            bp_r.append(pairlist[-1])

    # count all possible basepairs and all unpaired bases in the middle lists
    AU_m = bp_m.count(('A', 'U')) + bp_m.count(('U', 'A'))
    GU_m = bp_m.count(('G', 'U')) + bp_m.count(('U', 'G'))
    GC_m = bp_m.count(('G', 'C')) + bp_m.count(('C', 'G'))
    A_m = up_m.count(('A', '-'))
    C_m = up_m.count(('C', '-'))
    G_m = up_m.count(('G', '-'))
    U_m = up_m.count(('U', '-'))

    # normalize the individual values
    AU_m = float(AU_m) / len(bp_m)
    GU_m = float(GU_m) / len(bp_m)
    GC_m = float(GC_m) / len(bp_m)

    A_m = float(A_m) / len(up_m)
    C_m = float(C_m) / len(up_m)
    G_m = float(G_m) / len(up_m)
    U_m = float(U_m) / len(up_m)    
    
    results_m = [AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m]

    # count all possible basepairs and all unpaired bases in the edge lists 
    AU_r = bp_r.count(('A', 'U')) + bp_r.count(('U', 'A'))
    GU_r = bp_r.count(('G', 'U')) + bp_r.count(('U', 'G'))
    GC_r = bp_r.count(('G', 'C')) + bp_r.count(('C', 'G'))    

    A_r5 = up_r5.count(('A', '-'))
    C_r5 = up_r5.count(('C', '-'))
    G_r5 = up_r5.count(('G', '-'))
    U_r5 = up_r5.count(('U', '-'))

    A_r3 = up_r3.count(('A', '-'))
    C_r3 = up_r3.count(('C', '-'))
    G_r3 = up_r3.count(('G', '-'))
    U_r3 = up_r3.count(('U', '-'))

    # normalize the individual values
    AU_r = float(AU_r) / len(bp_r)
    GU_r = float(GU_r) / len(bp_r)
    GC_r = float(GC_r) / len(bp_r)

    A_r5 = float(A_r5) / len(up_r5)
    C_r5 = float(C_r5) / len(up_r5)
    G_r5 = float(G_r5) / len(up_r5)
    U_r5 = float(U_r5) / len(up_r5)    

    A_r3 = float(A_r3) / len(up_r3)
    C_r3 = float(C_r3) / len(up_r3)
    G_r3 = float(G_r3) / len(up_r3)
    U_r3 = float(U_r3) / len(up_r3)    
    
    results_r = [AU_r, GU_r, GC_r, A_r5, C_r5, G_r5, U_r5, A_r3, C_r3, G_r3, U_r3]
    
    results = results_m + results_r
    results.append(GC(sequence))
    return results

def main(): 
    """Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""
    print 'sequence, structure, energy_of_struct1, energy_of_struct2, AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m, AU_r, GU_r, GC_r, A_r5, C_r5, G_r5, U_r5, A_r3, C_r3, G_r3, U_r3, GCcontent, delta_energy'
    parser = argument_parser(sys.argv[1:])
    temperature_1, temperature_2 = parser.temperatures
    for line in sys.stdin:
        sequence, structure = line.split()
        delta_energy, energy_of_struct1, energy_of_struct2 = temperature_reactivity(sequence, structure, temperature_1, temperature_2)
        [AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m, AU_r, GU_r, GC_r, A_r5, C_r5, G_r5, U_r5, A_r3, C_r3, G_r3, U_r3, GCcontent] = bp_stats(sequence, structure)
        results = [sequence, structure, energy_of_struct1, energy_of_struct2, AU_m, GU_m, GC_m, A_m, C_m, G_m, U_m, AU_r, GU_r, GC_r, A_r5, C_r5, G_r5, U_r5, A_r3, C_r3, G_r3, U_r3, GCcontent, delta_energy]
        print ', '.join(map(str, results))
            
    
if __name__ == '__main__':
    sys.exit(main())

