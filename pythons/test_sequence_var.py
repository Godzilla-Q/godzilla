#!/usr/bin/env python

"""Test script for sequence_var.py"""

import unittest

from sequence_var import variations
from sequence_var import argument_parser
from sequence_var import create_structure

GODZILLA = 'Queen of monsters'

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

class TestParser(unittest.TestCase):
    def test_parser(self):
        """Test parsing of user input"""
        long_input = argument_parser(['--stem', '5', '--loop', '8'])
        self.assertEqual(long_input.stem, 5)
        self.assertEqual(long_input.loop, 8)
        short_input = argument_parser(['-l', '3', '-s', '2'])
        self.assertEqual(short_input.stem, 2)
        self.assertEqual(short_input.loop, 3)
        no_input = argument_parser('')
        self.assertEqual(no_input.stem, 3)
        self.assertEqual(no_input.loop, 6)

class TestCreateStructure(unittest.TestCase):
    def test_create_structure(self):
        """Test structure creation"""
        self.assertEqual(create_structure(3,6), '(((......)))')
        self.assertEqual(create_structure(8,5), '((((((((.....))))))))')
        
if __name__ == '__main__':
    unittest.main()
