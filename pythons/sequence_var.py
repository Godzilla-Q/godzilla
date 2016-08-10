#!/usr/bin/env python

"""Generate new sequence variation: switches the first base pair of a structure"""

# from __future__ import print_function
import sys
import RNA
import re
from Bio import SeqIO
from Bio.Seq import Seq
from ss_dotplot import versions_used

"""Switches the first base pair of the first stem loop of a sequence"""

def determine_bp( structure ):
    # gives out the position of the first basepair in a structure
    ptable = RNA.ptable(structure)[1:]
    a, b = 0, 0
    for index, item in enumerate(ptable):
        if item:
            a = index + 1
            b = item
            break
    return a, b

def second__bp( structure ):
    # gives out the position of the first basepair in a structure
    ptable = RNA.ptable(structure)[1:]
    a, b = 0, 0
    for index, item in enumerate(ptable):
        if item:
            a = index + 2
            b = item - 1
            break
    return a, b

def switch_bp( seqrec, a, b ):
    # switches basepair in SeqRecord indicated by the positions a and b to their complement
    complement = seqrec.reverse_complement()[::-1]
    if a < b:
        new = seqrec[:a-1] + complement[a-1] + seqrec[a:b-1] + complement[b-1]+ seqrec[b:]
    else:
        new = seqrec
        print 'WARNING: No valid basepairs in this structure'
    return new

def oldmain():
    for original in SeqIO.parse(sys.stdin, "fasta"):
        # obtain the sequence string
        sequ = str(original.seq)
        # generate the fold compound
        fc = RNA.fold_compound(sequ)
        # calculate mfe structure and mfe
        struct, mfe = fc.mfe()
        # change base pair
        a, b = determine_bp(struct)
        new = switch_bp(original,a,b)
        print new.seq
    print versions_used()

def calculate_values( sequence, structure, temperature1, temperature2 ):
    """Evaluate energy, entropy, enthalpy for each temperature"""
    temperature_model1 = RNA.md()
    temperature_model2 = RNA.md()
    temperature_model1.temperature = temperature1
    temperature_model2.temperature = temperature2
    fc1 = RNA.fold_compound(sequence, temperature_model1)
    fc2 = RNA.fold_compound(sequence, temperature_model2)
    energy_of_struct1 = fc1.eval_structure(structure)
    energy_of_struct2 = fc2.eval_structure(structure)
    delta_energy = energy_of_struct2 - energy_of_struct1
    delta_entropy = (-delta_energy/(temperature2 - temperature1))
    delta_enthalpy = (energy_of_struct1 + temperature1 * delta_entropy)
#    delta_enthalpy2 = (energy_of_struct2 + temperature2 * delta_entropy)    
    return (delta_energy, delta_entropy, delta_enthalpy)
    
def main():
    seq_object = Seq('ccccaaaaaagggg')
    structure = '((((......))))'
    results = calculate_values(str(seq_object), structure, 30, 37)
    a, b = determine_bp(structure)
    new_object = switch_bp(seq_object,a,b)
    results_new = calculate_values(str(new_object), structure, 30, 37)
    print str(seq_object), results
    print str(new_object), results_new

if __name__ == '__main__':
    sys.exit(main())
 

    

