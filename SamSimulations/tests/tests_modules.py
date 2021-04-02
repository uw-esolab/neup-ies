#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

import unittest, os,sys
from util.FileMethods import FileMethods

class TestPySAMModules(unittest.TestCase):
    """
    Unit tests for PySAM module setup
    """
    module_names = ['tcsmolten_salt',
                         'nuclear_tes']
    
    module_tests = ['generic_mspt',
                         'generic_nuclear_model1']
    
    
    def setUp(self):
        return 0
        