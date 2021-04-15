#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

from modules.NuclearTES import NuclearTES
import unittest

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
        
        # creating instance of module
        self.nuctes = NuclearTES(json_name='tests/test_nuctes')
        self.nuctes_name = self.nuctes.__class__.__name__
        
    
    def tearDown(self):
        """ Deleting instance of NuclearTES at end of each test
        """
        
        # deleting each module
        del self.nuctes
        del self.nuctes_name
    

    def test_store_csv_arrays(self):
        """ Testing the storage of csv arrays in NuclearTES
        
        NOTE: this assumes that the method was already called in __init__
        """
        
        # in __init__ we call self.store_csv_arrays which is overloaded in NucTes
        #    double checking that the NuclearTES-specific attributes are generated
        attr_list = ['df_array' , 'ud_array', 'wl_array', 'hp_array', 
                     'gc_array' , 'em_array', 'fm_array' ]
    
        
        # looping through all defined attributes
        for attr in attr_list:
            #TODO: check that the filepaths exist from PySAM_dict?
            
            # checking that the attribute exists
            self.assertTrue(hasattr(self.nuctes,attr) ,
                            "Something went wrong when {0} called 'store_csv_arrays method, {1} does not exist".format(self.nuctes_name, attr) )
            
if __name__ == "__main__":
    unittest.main()
            