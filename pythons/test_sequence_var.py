#!/usr/bin/env python

"""Test script for sequence_var.py"""

import unittest
from sequence_var import temperature_reactivity

import sys
import RNA
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from ss_dotplot import versions_used
from itertools import product
from operator import itemgetter

#auaaaaaaauau (((......))) 30 37
#0.430000066757 -0.0614285809653 0.927142551967

class TestReactivityValues(unittest.TestCase):

    def test_temperature_reactivity(self):
        self.assertEqual(temperature_reactivity('auaaaaaaauau', '(((......)))', 30, 37), (0.43000006675720215, -0.061428580965314596, 0.9271425519670757))
        self.assertEqual(temperature_reactivity('aaaaaaaaaauuuu', '((((......))))', 30, 37), (0.5, -0.07142857142857142, 0.25714295251028885))
        self.assertEqual(temperature_reactivity('gggaaaaaauuu', '(((......)))', 30, 37), (0.6999998092651367, -0.0999999727521624, 1.4000009128025597))
        self.assertEqual(temperature_reactivity('gcccaaaaaagggc', '((((......))))', 30, 37), (1.0500001907348633, -0.15000002724783762, -11.650000912802561))

if __name__ == '__main__':
    unittest.main()
