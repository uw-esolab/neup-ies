#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 14:19:23 2021

@author: gabrielsoto
"""

import unittest, os
import numpy as np
from util.FileMethods import FileMethods

class TestUtils(unittest.TestCase):
    """
    Unit tests for PySAM utils
    """
    
    def test_csvread(self):
        
        #testing with test1.csv file, only has 1 column
        csvpath1 = os.path.join( FileMethods.neup_dir , "data/tests/test1.csv")
        data_array1 = FileMethods.read_csv_through_pandas(csvpath1)
        
        #----checking that the output array is 1 dimensional as expected
        self.assertEqual(type(data_array1),np.ndarray, "test1 output is not an ndarray")
        self.assertEqual(len(data_array1.shape),1,     "test1 is not a 1D array") 
        
        #testing with test1.csv file, has 2 columns
        csvpath2 = os.path.join( FileMethods.neup_dir , "data/tests/test2.csv")
        data_array2 = FileMethods.read_csv_through_pandas(csvpath2)
        
        #----checking that the output array is a list of lists
        self.assertTrue(isinstance(data_array2,        list),  "test2 output is not a list" )
        self.assertTrue(isinstance(data_array2[0],     list),  "test2 output is not a list of lists" )
        self.assertFalse(isinstance(data_array2[0][0], list),  "test2 output has an extra column?" )
