#!/usr/bin/env python

"""Test script for sequence_var.py"""

import unittest
from sequence_var import temperature_reactivity
from sequence_var import variations
from sequence_var import sort_results
from sequence_var import sequence_properties
from sequence_var import argument_parser
from sequence_var import create_structure
from sequence_var import energy_profile

import sys
import RNA
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from ss_dotplot import versions_used
from itertools import product
from operator import itemgetter

GODZILLA = 'Queen of monsters'

class TestReactivityValues(unittest.TestCase):

    def test_temperature_reactivity(self):
        "Compare to previously calculated results"
        self.assertEqual(temperature_reactivity('auaaaaaaauau', '(((......)))', 30, 37), (0.43000006675720215, -0.061428580965314596, 0.9271425519670757))
        self.assertEqual(temperature_reactivity('aaaaaaaaaauuuu', '((((......))))', 30, 37), (0.5, -0.07142857142857142, 0.25714295251028885))
        self.assertEqual(temperature_reactivity('gggaaaaaauuu', '(((......)))', 30, 37), (0.6999998092651367, -0.0999999727521624, 1.4000009128025597))
        self.assertEqual(temperature_reactivity('gcccaaaaaagggc', '((((......))))', 30, 37), (1.0500001907348633, -0.15000002724783762, -11.650000912802561))

class TestVariations(unittest.TestCase):
    example_variations = ['aaaaauu',
                          'auaaaau',
                          'agaaacu',
                          'acaaagu',
                          'agaaauu',
                          'auaaagu',
                          'uaaaaua',
                          'uuaaaaa',
                          'ugaaaca',
                          'ucaaaga',
                          'ugaaaua',
                          'uuaaaga',
                          'gaaaauc',
                          'guaaaac',
                          'ggaaacc',
                          'gcaaagc',
                          'ggaaauc',
                          'guaaagc',
                          'caaaaug',
                          'cuaaaag',
                          'cgaaacg',
                          'ccaaagg',
                          'cgaaaug',
                          'cuaaagg',
                          'gaaaauu',
                          'guaaaau',
                          'ggaaacu',
                          'gcaaagu',
                          'ggaaauu',
                          'guaaagu',
                          'uaaaaug',
                          'uuaaaag',
                          'ugaaacg',
                          'ucaaagg',
                          'ugaaaug',
                          'uuaaagg']

    def test_variations(self):
        "Check number of variations"
        self.assertEqual(len(list(variations(3,1))), 216)
        self.assertEqual(len(list(variations(4,20))), 1296)
        self.assertEqual(list(variations(2,3)), self.example_variations)

class TestSortResults(unittest.TestCase):
    list_of_tuples_sortby_0 = [(0, 4, 3, 2, 1),
                               (1, 0, 4, 3, 2),
                               (2, 1, 0, 4, 3),
                               (3, 2, 1, 0, 4),
                               (4, 3, 2, 1, 0)]

    list_of_tuples_sortby_1 = [(1, 0, 4, 3, 2),
                               (2, 1, 0, 4, 3),
                               (3, 2, 1, 0, 4),
                               (4, 3, 2, 1, 0),
                               (0, 4, 3, 2, 1)]

    list_of_tuples_sortby_2 = [(2, 1, 0, 4, 3),
                               (3, 2, 1, 0, 4),
                               (4, 3, 2, 1, 0),
                               (0, 4, 3, 2, 1),
                               (1, 0, 4, 3, 2)]

    list_of_tuples_sortby_3 = [(3, 2, 1, 0, 4),
                               (4, 3, 2, 1, 0),
                               (0, 4, 3, 2, 1),
                               (1, 0, 4, 3, 2),
                               (2, 1, 0, 4, 3)]

    list_of_tuples_sortby_4 = [(4, 3, 2, 1, 0),
                               (0, 4, 3, 2, 1),
                               (1, 0, 4, 3, 2),
                               (2, 1, 0, 4, 3),
                               (3, 2, 1, 0, 4)]     
    
    def test_sort_results(self):
        """Check sorting"""
        self.assertEqual(sort_results(self.list_of_tuples_sortby_0, 0), self.list_of_tuples_sortby_0)
        self.assertEqual(sort_results(self.list_of_tuples_sortby_0, 1), self.list_of_tuples_sortby_1)
        self.assertEqual(sort_results(self.list_of_tuples_sortby_0, 2), self.list_of_tuples_sortby_2)
        self.assertEqual(sort_results(self.list_of_tuples_sortby_0, 3), self.list_of_tuples_sortby_3)
        self.assertEqual(sort_results(self.list_of_tuples_sortby_0, 4), self.list_of_tuples_sortby_4)

class TestSequenceProperties(unittest.TestCase):
    def test_sequence_properties(self):
        """Check sequence properties results"""
        self.assertTrue(sequence_properties('ccggg') == 100)
        self.assertTrue(sequence_properties('ccgga') == 80)
        self.assertTrue(sequence_properties('ccgau') == 60)
        self.assertTrue(sequence_properties('ccaua') == 40)
        self.assertTrue(sequence_properties('cuaua') == 20)
        self.assertTrue(sequence_properties('auaua') == 0)
        

class TestParser(unittest.TestCase):
    def test_parser(self):
        """Test parsing of user input"""
        long_input = argument_parser(['--stem', '5', '--temperatures', '20', '40', '--loop', '8'])
        self.assertEqual(long_input.stem, 5)
        self.assertEqual(long_input.loop, 8)
        self.assertEqual(long_input.temperatures, [20, 40])
        short_input = argument_parser(['-l', '3', '-s', '2', '-T', '12', '80'])
        self.assertEqual(short_input.stem, 2)
        self.assertEqual(short_input.loop, 3)
        self.assertEqual(short_input.temperatures, [12, 80])
        no_input = argument_parser('')
        self.assertEqual(no_input.stem, 3)
        self.assertEqual(no_input.loop, 6)
        self.assertEqual(no_input.temperatures, [30, 37])

class TestCreateStructure(unittest.TestCase):
    def test_create_structure(self):
        """Test structure creation"""
        self.assertEqual(create_structure(3,6), '(((......)))')
        self.assertEqual(create_structure(8,5), '((((((((.....))))))))')

class TestEnergyProfile(unittest.TestCase):
    example_sequence = 'cgcaaagcg'
    example_structure = '(((...)))'
    example_result = (0, -2.4, -3.4, 5.4)
    
    def test_energy_profile(self):
        """Test evaluation of base pair energy"""
        # maybe include more tests for strange cases?
        self.assertEqual(energy_profile(self.example_sequence, self.example_structure), self.example_result)
        
if __name__ == '__main__':
    unittest.main()
