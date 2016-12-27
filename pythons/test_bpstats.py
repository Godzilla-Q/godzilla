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

        expected_results1 = [0.0025,
                             0.004375,
                             0.003125,
                             0.002903225806451613,
                             0.0035483870967741938,
                             0.001935483870967742,
                             0.0016129032258064516,
                             0.006,
                             0.0,
                             0.004,
                             0.004117647058823529,
                             0.000588235294117647,
                             0.0029411764705882353,
                             0.002352941176470588,
                             48.0]
        self.assertEqual(bp_stats(sequence1, structure1), expected_results1)

        structure2 = '.(((......))).'
        sequence2 = 'AGUGCUCCCACACC'

        expected_results2 = [0.07142857142857142,
                             0.0,
                             0.0,
                             0.0,
                             0.05357142857142857,
                             0.0,
                             0.017857142857142856,
                             0.0,
                             0.0,
                             0.07142857142857142,
                             0.03571428571428571,
                             0.03571428571428571,
                             0.0,
                             0.0,
                             64.28571428571429]

        self.assertEqual(bp_stats(sequence2, structure2), expected_results2)

        structure3 = '(((......)))'
        sequence3 = 'CCCAAAAAAGGG'

        expected_results3 = [0.0, 0.0, 0.08333333333333333, 0.08333333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08333333333333333, 0.08333333333333333, 0.0, 0.0, 0.0, 50.0]

        self.assertEqual(bp_stats(sequence3, structure3), expected_results3)
        

    

        
if __name__ == '__main__':
    unittest.main()
