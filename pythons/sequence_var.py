#!/usr/bin/env python

"""Generate variations of the bases of the stem; calculate the difference in thermodynamic parameters for two different temperatures"""

import sys
import RNA
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from ss_dotplot import versions_used
from itertools import product
from operator import itemgetter
from subprocess import Popen, PIPE


def temperature_reactivity( sequence, structure, temperature1, temperature2 ):
    """Evaluate temperature-dependent difference in energy, entropy, enthalpy."""
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

def variations( stem_length, loop_length ):
    """Yield all sequence variations for a stem loop with fixed loop"""
    basepair_set = ['au', 'ua', 'gc', 'cg', 'gu', 'ug']
    center = 'a' * loop_length    
    for basepairs in product(basepair_set, repeat = stem_length):
        left, right = zip(*basepairs)
        sequence = (''.join(left), center, ''.join(right[::-1])) 
        yield ''.join(sequence)

def sort_results( list_of_tuples, number ):
    """Sort list of tuples with respect to the attribute number 'number'"""
    return sorted(list_of_tuples, key=itemgetter(number))

def sequence_properties( sequence ):
    """Calculate GC content"""
    return GC(sequence)

def argument_parser( arguments ):
    """If arguments from command line are present, override default parameters. Ensure that sequence matches structure length. Manage help messages"""
    
    parser = argparse.ArgumentParser(
        description='Generate variations of the bases of the stem; calculate the difference in thermodynamic parameters for two different temperatures.',
        epilog='All hail Godzilla, Queen of monsters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-T', '--temperatures',
        nargs=2,
        default=[30, 37],
        type=int,
        help='two temperatures in degree Celsius',
        metavar=('TEMPERATURE1', 'TEMPERATURE2')
    )
    parser.add_argument(
        '-s', '--stem',
         default=3,
        type=int,
        choices=xrange(2, 6),
        help='length of stem'
    )
    parser.add_argument(
        '-l', '--loop',
         default=6,
        type=int,
        choices=xrange(3, 10),
        help='length of loop'
    )

    return parser.parse_args(arguments)

def create_structure( stem, loop ):
    """Create structure string corresponding to the sequences."""
    return ''.join(['('*stem, '.'*loop, ')'*stem])

def energy_profile( sequence, structure ):
    """Calculate energy contribution of each basepair - call RNAeval from the shell and to obtain the values from the stdout of RNAeval manually"""
    eval_input = '\n'.join([sequence,structure])
    subprocess_eval = Popen(['RNAeval', '-v'], stdin=PIPE, stdout=PIPE, universal_newlines=True)
    eval_output = subprocess_eval.communicate(input=eval_input)
    result_list = []
    for line in eval_output[0].splitlines():
        if ':' in line:
            trash, treasure = line.split(':')
            result_list.append(float(treasure)/100)
    return tuple(result_list)
    # check the units of the results! 10cal/mol? Any way to improve my code/make it more elegant?

def main():
    """Generate variations of the bases of the stem; calculate the difference in
thermodynamic parameters for two different temperatures."""
    parser = argument_parser(sys.argv[1:])
    temperature_1, temperature_2 = parser.temperatures
    stem = parser.stem
    loop = parser.loop
    structure = create_structure(stem, loop)
    results = []
    for sequence in variations(stem, loop):
        energy, entropy, enthalpy = temperature_reactivity(sequence, structure, temperature_1, temperature_2)
        gc_content = sequence_properties(sequence)
        e_profile = energy_profile(sequence, structure)
        results.append((sequence, energy, entropy, enthalpy, gc_content, e_profile))
    for number in range(1, 4):
        for line in sort_results(results, number):
            print line[0], line[1], line[2], line[3], line[4], line[5]
        print '\n'

    
if __name__ == '__main__':
    sys.exit(main())
