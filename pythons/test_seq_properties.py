#!/usr/bin/env python

"""Test script for seq_properties.py"""

import unittest

from seq_properties import temperature_reactivity
from seq_properties import gc_content
from seq_properties import argument_parser
from seq_properties import energy_profile

GODZILLA = 'Queen of monsters'

class TestReactivityValues(unittest.TestCase):
    def test_temperature_reactivity(self):
        "Compare to previously calculated results"
        self.assertEqual(temperature_reactivity('auaaaaaaauau',
                                                '(((......)))',
                                                30,
                                                37), (
                                                    0.13437501885928185,
                                                    2.7699999809265137,
                                                    3.200000047683716))
        self.assertEqual(temperature_reactivity('aaaaaaaaaauuuu',
                                                '((((......))))',
                                                30,
                                                37), (
                                                    0.17241378743356547,
                                                    2.4000000953674316,
                                                    2.9000000953674316))
        self.assertEqual(temperature_reactivity('gggaaaaaauuu',
                                                '(((......)))',
                                                30,
                                                37), (
                                                    0.1372548671283884,
                                                    4.400000095367432,
                                                    5.099999904632568))
        self.assertEqual(temperature_reactivity('gcccaaaaaagggc',
                                                '((((......))))',
                                                30,
                                                37), (
                                                    0.17213118150009377,
                                                    -7.150000095367432,
                                                    -6.099999904632568))

class TestGCContent(unittest.TestCase):
    def test_gc_content(self):
        """Check sequence properties results"""
        self.assertTrue(gc_content('ccggg') == 100)
        self.assertTrue(gc_content('ccgga') == 80)
        self.assertTrue(gc_content('ccgau') == 60)
        self.assertTrue(gc_content('ccaua') == 40)
        self.assertTrue(gc_content('cuaua') == 20)
        self.assertTrue(gc_content('auaua') == 0)
        
class TestParser(unittest.TestCase):
    def test_parser(self):
        """Test parsing of user input"""
        long_input = argument_parser(['--temperatures', '20', '40'])
        self.assertEqual(long_input.temperatures, [20, 40])
        short_input = argument_parser(['-T', '12', '80'])
        self.assertEqual(short_input.temperatures, [12, 80])
        no_input = argument_parser('')
        self.assertEqual(no_input.temperatures, [30, 37])

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
