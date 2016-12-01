#!/usr/bin/env python

"""Test script for sequence_var.py"""

import unittest

from random_sequences_pipeline import argument_parser
from random_sequences_pipeline import unconstrained_structure
from random_sequences_pipeline import random_sequences
from random_sequences_pipeline import fold_sequences

GODZILLA = 'Queen of monsters'

class TestParser(unittest.TestCase):
    def test_parser(self):
        """Test parsing of user input"""
        long_input = argument_parser(['--length', '20',
                                      '--number1', '9',
                                      '--seed1', '100',
                                      '--number2', '8',
                                      '--seed2', '99'])
        self.assertEqual(long_input.length, 20)
        self.assertEqual(long_input.number1, 9)
        self.assertEqual(long_input.seed1, 100)
        self.assertEqual(long_input.number2, 8)
        self.assertEqual(long_input.seed2, 99)       
        short_input = argument_parser(['-l', '90',
                                      '-n1', '5',
                                      '-s1', '80',
                                      '-n2', '3',
                                      '-s2', '1000'])
        self.assertEqual(short_input.length, 90)
        self.assertEqual(short_input.number1, 5)
        self.assertEqual(short_input.seed1, 80)
        self.assertEqual(short_input.number2, 3)
        self.assertEqual(short_input.seed2, 1000)       
        no_input = argument_parser('')
        self.assertEqual(no_input.length, 100)
        self.assertEqual(no_input.number1, 10)
        self.assertEqual(no_input.seed1, 1)
        self.assertEqual(no_input.number2, 100)
        self.assertEqual(no_input.seed2, 1)       

class TestUnconstrainedStruct(unittest.TestCase):
    def test_unconstrained_struct(self):
        self.assertEqual(unconstrained_structure(28), '............................')
        self.assertEqual(unconstrained_structure(0), '')
        self.assertEqual(unconstrained_structure(1), '.')
        
class TestRandomSeq(unittest.TestCase):
    def test_random_seq(self):
        """Test structure creation"""
        expected_result1 = ['CCUACAGAGUGGUCAUUAGGAAAGUCAC',
                           'UUUGGUGACUUCCAGCUUACAUAGGGCC',
                           'AAGCGAGCUUAACAUUCAUGCGCCGUGA']
        actual_result1 = list(random_sequences('............................', 3, 5))
        self.assertEqual(expected_result1, actual_result1)
        expected_result2 = ['CUGGAUCGGUGUUUAGAGUAGUCACAUU',
                            'GCUGAACCUAUUAGGUGGAGCCAACUCU',
                            'GUCUUUAGCGAUGGAUCGUUAAUAGAUG',
                            'CUGCGACCCACAUCGGGUCUAAUUGGAU']
        actual_result2 = list(random_sequences('(((..........)))(((......)))', 4, 5))        
        self.assertEqual(expected_result2, actual_result2)
        
class TestRnaFold(unittest.TestCase):
    def test_rna_fold(self):
        sequence = 'CCUACAGAGUGGUCAUUAGGAAAGUCACUUUGGUGACUUCCAGCUUACAUAGGGCAAGCGAGCUUAACAUU'
        expected_result = '.((((...))))......((((.(((((....)))))))))(((((.(.........).))))).......'
        actual_result = fold_sequences(sequence)
        self.assertEqual(expected_result, actual_result)
        
if __name__ == '__main__':
    unittest.main()
