#!/usr/bin/env python

"""Generate a number of random sequences, fold them, and generate a number of random sequences for each structure"""

import sys
import argparse
import RNA
import RNAblueprint as rd
from itertools import product

def argument_parser( arguments ):
    """If arguments from command line are present, override default parameters. Ensure that sequence matches structure length. Manage help messages"""
    
    parser = argparse.ArgumentParser(
        description='Generate a number of random sequences, fold them, and generate a number of random sequences for each structure.',
        epilog='All hail Godzilla, Queen of monsters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-l', '--length',
         default=100,
        type=int,
        choices=xrange(1, 100),
        help='length of random sequences'
    )
    parser.add_argument(
        '-n1', '--number1',
         default=10,
        type=int,
        choices=xrange(1, 10),
        help='number of random sequences created at the start'
    )
    parser.add_argument(
        '-s1', '--seed1',
         default=1,
        type=int,
        help='random number generator seed for random sequences'
    )
    parser.add_argument(
        '-n2', '--number2',
         default=100,
        type=int,
        choices=xrange(1, 100),
        help='number of random sequences for each structure'
    )
    parser.add_argument(
        '-s2', '--seed2',
         default=1,
        type=int,
        help='random number generator seed for sequences for each structure'
    )

    return parser.parse_args(arguments)

def unconstrained_structure( length ):
    """Create structure string corresponding to the sequences."""
    return '.' * length

def random_sequences( structure, number, seed ):
    """Create random sequences compatible with a given structure."""
    # construct dependency graph with RNAblueprint
    dg = rd.DependencyGraphMT([structure], '', seed)
    # create sequence variations and print them
    for i in range (0, number):
        dg.sample()
        yield dg.get_sequence()

def fold_sequences( sequence ):
    """Fold random sequence"""
    fc = RNA.fold_compound(sequence)
    (structure, mfe) = fc.mfe()
    return structure

def main():
    """Generate a stem-loop structure and all possible variations of the bases of the stem."""
    parser = argument_parser(sys.argv[1:])
    length = parser.length
    n1 = parser.n1
    s1 = parser.s1
    n2 = parser.n2
    s2 = parser.s2
    structure = unconstrained_structure(length)
    random_structures = []
    for randseq in random_sequences(structure, n1, s1):
        random_structures.append(fold_sequences(randseq))
    for struct in random_structures:
        for seq in random_sequences(struct, n2, s2):
            print seq, struct
    
if __name__ == '__main__':
    sys.exit(main())
