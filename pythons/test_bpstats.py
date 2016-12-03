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

        expected_results1 = [0.08, 0.14, 0.10, 0.09, 0.11, 0.06, 0.05, 0.12, 0.0, 0.08, 0.07, 0.01, 0.05, 0.04, 48.0]
        self.assertEqual(bp_stats(sequence1, structure1), expected_results1)

        structure2 = '.(((......))).'
        sequence2 = 'AGUGCUCCCACACC'

        expected_results2 = [0.14285714285714285, 0.0, 0.0, 0.0, 0.21428571428571427, 0.0, 0.07142857142857142, 0.0, 0.0, 0.2857142857142857, 0.14285714285714285, 0.14285714285714285, 0.0, 0.0, 64.28571428571429]

        self.assertEqual(bp_stats(sequence2, structure2), expected_results2)
        

    

        
if __name__ == '__main__':
    unittest.main()
