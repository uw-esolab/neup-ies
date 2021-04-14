#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

from modules.NuclearTES import NuclearTES
import unittest, os

class TestNuclearTES(unittest.TestCase):
    """
    Unit tests for Nuclear TES module
    
    This testing suite is meant to test the nuclear TES module meant to
    represent Model 1 of the NE-2 project. This models an LFR plant supplying 
    power to an sCO2 power block with thermal energy storage (TES) in the form 
    of molten salt tanks. 
    """

    def setUp(self):
        """ Creating instance of NuclearTES upon start of each test
        """
        
        self.nuctes = NuclearTES(json_name='tests/test_nuctes')
        
    
    def tearDown(self):
        """ Deleting instance of NuclearTES at end of each test
        """
        
        # deleting each module
        del self.nuctes