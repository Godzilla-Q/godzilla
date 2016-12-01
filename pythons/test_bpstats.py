#!/usr/bin/env python

"""Test script for seq_properties.py"""

import unittest

from bpstats import bp_stats

GODZILLA = 'Queen of monsters'

class TestBPStats(unittest.TestCase):
    def test_bp_stats(self):
        "Compare to previously calculated results"

        self.assertTrue(bp_stats('AAAAA', '.....') == [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        self.assertTrue(bp_stats('CCCCC', '.....') == [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 100.0])
        self.assertTrue(bp_stats('GGGGG', '.....') == [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 100.0])
        self.assertTrue(bp_stats('UUUUU', '.....') == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
        self.assertTrue(bp_stats('AUAGGGGGGUAU', '(((......)))') == [0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 50.0])
        self.assertTrue(bp_stats('GGGAAAAAAUUU', '(((......)))') == [0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 25.0])
        self.assertTrue(bp_stats('GGGAAAAAACCC', '(((......)))') == [0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 50.0])
        self.assertFalse(bp_stats('AAAAA', '.....') == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        
if __name__ == '__main__':
    unittest.main()
