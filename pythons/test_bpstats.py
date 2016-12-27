#!/usr/bin/env python

"""Test script for seq_properties.py"""

import unittest

from bpstats import bp_stats_old
from bpstats import bp_stats

GODZILLA = 'Queen of monsters'

class TestBPStats(unittest.TestCase):
    def test_bp_stats_old(self):
        "Compare to previously calculated results"

        self.assertEqual(bp_stats_old('AAAAA', '.....'), [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.assertEqual(bp_stats_old('CCCCC', '.....'), [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 100.0])
        self.assertEqual(bp_stats_old('GGGGG', '.....'), [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 100.0])
        self.assertEqual(bp_stats_old('UUUUU', '.....'), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
        self.assertEqual(bp_stats_old('AUAGGGGGGUAU', '(((......)))'), [0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 50.0])
        self.assertEqual(bp_stats_old('GGGAAAAAAUUU', '(((......)))'), [0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 25.0])
        self.assertEqual(bp_stats_old('GGGAAAAAACCC', '(((......)))'), [0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 50.0])
        self.assertFalse(bp_stats_old('AAAAA', '.....') == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    def test_bp_stats(self):
        structure1 = '....((((.....))))........((((...))))(((((((.(((((((.......)))))))...)))))))....((((.......))))......'
        sequence1 = 'AAGUCGAAACCGGUUUGUACUGCCGGGUCAAAGACCUGUCGUUAAAGUUUUAUCGCGGAGGGUUUUAAAGCGGCAAACGCGAUUUAGUCGAUCGCCUCAA'

        expected_results1 = [0.25,
                             0.4375,
                             0.3125,
                             0.2903225806451613,
                             0.35483870967741938,
                             0.1935483870967742,
                             0.16129032258064516,
                             0.6,
                             0.0,
                             0.4,
                             0.3333333333333333,
                             0.0,
                             0.5555555555555556,
                             0.1111111111111111,
                             0.5555555555555556,
                             0.1111111111111111,
                             0.0,
                             0.3333333333333333,
                             48.0]
        
        self.assertEqual(bp_stats(sequence1, structure1), expected_results1)

        structure2 = '.(((......))).'
        sequence2 = 'AGUGCUCCCACACC'

        expected_results2 = [1.0,
                             0.0,
                             0.0,
                             0.0,
                             0.75,
                             0.0,
                             0.25,
                             0.0,
                             0.0,
                             1.0,
                             1.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             1.0,
                             0.0,
                             0.0,
                             64.28571428571429]

        self.assertEqual(bp_stats(sequence2, structure2), expected_results2)

        structure3 = '(((......)))'
        sequence3 = 'CCCAAAAAAGGG'

        expected_results3 = [0.0,
                             0.0,
                             1.0,
                             1.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             0.0,
                             1.0,
                             1.0,
                             0.0,
                             0.0,
                             0.0,
                             1.0,
                             0.0,
                             0.0,
                             0.0,
                             50.0]

        self.assertEqual(bp_stats(sequence3, structure3), expected_results3)
        

    

        
if __name__ == '__main__':
    unittest.main()
