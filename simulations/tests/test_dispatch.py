#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:54:26 2021

@author: gabrielsoto
"""


from dispatch.GeneralDispatch import GeneralDispatch
import unittest

class TestGeneralDispatch(unittest.TestCase):
    """
    Unit tests for General Dispatch module
    
    This testing suite is meant to test the general Dispatch module of 
    the NE2 architecture. Here we have 
    """

    def setUp(self):
        """ Creating instance of GeneralDispatch upon start of each test
        """
        
        # creating instance of module
        self.dispatch = GeneralDispatch({})
        
    
    def tearDown(self):
        """ Deleting instance of GeneralDispatch at end of each test
        """
        
        # deleting each module
        del self.dispatch
    

if __name__ == "__main__":
    unittest.main()
            