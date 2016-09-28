#!/usr/bin/env python

"""Generate a stem-loop structure and all possible variations of the bases of the stem"""

import sys
import argparse
from itertools import product

def variations( stem_length, loop_length ):
    """Yield all sequence variations for a stem loop with fixed loop"""
    basepair_set = ['au', 'ua', 'gc', 'cg', 'gu', 'ug']
    center = 'a' * loop_length    
    for basepairs in product(basepair_set, repeat = stem_length):
        left, right = zip(*basepairs)
        sequence = (''.join(left), center, ''.join(right[::-1])) 
        yield ''.join(sequence)

def argument_parser( arguments ):
    """If arguments from command line are present, override default parameters. Ensure that sequence matches structure length. Manage help messages"""
    
    parser = argparse.ArgumentParser(
        description='Generate a stem-loop structure and all possible variations of the bases of the stem.',
        epilog='All hail Godzilla, Queen of monsters.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
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

def main():
    """Generate a stem-loop structure and all possible variations of the bases of the stem."""
    parser = argument_parser(sys.argv[1:])
    stem = parser.stem
    loop = parser.loop
    structure = create_structure(stem, loop)
    results = []
    for sequence in variations(stem, loop):
        print sequence, structure
    
if __name__ == '__main__':
    sys.exit(main())
