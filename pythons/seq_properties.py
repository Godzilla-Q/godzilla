#!/usr/bin/env python

"""Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""

import sys
import RNA
import argparse
from Bio.SeqUtils import GC
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
    # normalize delta_energy, add 0.001 to prevent division by 0
    delta_energy = abs((energy_of_struct2 - energy_of_struct1) / (energy_of_struct2 + 0.001))

    return (delta_energy, energy_of_struct1, energy_of_struct2)

def gc_content( sequence ):
    """Calculate GC content"""
    return GC(sequence)

def argument_parser( arguments ):
    """If arguments from command line are present, override default parameters. Ensure that sequence matches structure length. Manage help messages"""
    
    parser = argparse.ArgumentParser(
        description='Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence',
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
    return parser.parse_args(arguments)

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

def main():
    """Calculate the difference in thermodynamic parameters for two different temperatures and calculate descriptors for the sequence"""
    print 'sequence, structure, energy_of_struct1, energy_of_struct2, gc_content, e0, e1, e2, e3, e4, e5, delta_energy'    
    parser = argument_parser(sys.argv[1:])
    temperature_1, temperature_2 = parser.temperatures
    for line in sys.stdin:
        sequence, structure = line.split()
        delta_energy, energy_of_struct1, energy_of_struct2 = temperature_reactivity(sequence, structure, temperature_1, temperature_2)
        gc = gc_content(sequence)
        e0, e1, e2, e3, e4, e5 = energy_profile(sequence, structure) #improve flexibility!
        print '{},{},{},{},{},{},{},{},{},{},{},{}'.format(sequence, structure, energy_of_struct1, energy_of_struct2, gc, e0, e1, e2, e3, e4, e5, delta_energy)

    
if __name__ == '__main__':
    sys.exit(main())
