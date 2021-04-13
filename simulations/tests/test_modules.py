#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

from modules.GenericSSCModule import GenericSSCModule
from modules.NuclearTES import NuclearTES
import unittest, os


class TestPySAMModules(unittest.TestCase):
    """
    Unit tests for PySAM module setup
    
    This testing suite is meant to test the parent class of all modules,
    that being GenericSSCModule. All child classes will share common methods
    and attributes which will be tested here. 
    
    For each individual class, there will be bespoke test classes.
    """

    # test that a full, single run is close to results from chainlinking
    
    # test that outputs are communicated from Plant to Grid to SO
    
    def setUp(self):
        """ Creating instances of modules upon start of each test
        """
        self.genmod = GenericSSCModule()
        self.nuctes = NuclearTES()
    
    
    def tearDown(self):
        """ Deleting instances of modules at end of each test
        """
        del self.genmod
        del self.nuctes


    def test__init__(self):
        """ Testing the shared processes in __init__ of all modules
        """
        
        # defining modules after calling setUp( ) method
        mod_list  = [self.genmod, self.nuctes]
        
        # defining necessary attributes to be declared in __init__
        attr_list = ['json_name' , 'plant_name', 'ssc_horizon', 'pyomo_horizon', 
                     'SSC_dict']
        
        # looping through all defined modules + attributes
        for mod in mod_list:
            for attr in attr_list:
                
                # checking that all attributes exist in every module
                self.assertTrue(hasattr(mod,attr) , 
                                "{0} does not have '{1}' attribute".format(mod.__class__.__name__,attr))
            
            # checking that the self.store_csv_arrays( ) method was called correctly
            self.assertTrue(hasattr(mod,"solar_resource_file") ,
                            "Something went wrong when {0} called 'store_csv_arrays' method".format(mod.__class__.__name__) )


    def test_store_csv_arrays(self):
        """ Testing the shared processes in __init__ of all modules
        """
        
        # defining modules after calling setUp( ) method
        mod_list  = [self.genmod, self.nuctes]
        
        # looping through all modules
        for mod in mod_list:
            
            # checking that the solar_resource_file path attribute exists
            self.assertTrue(hasattr(mod,"solar_resource_file") ,
                            "Something went wrong when {0} called 'store_csv_arrays' method".format(mod.__class__.__name__) )
            
            # checking that the actual filepath exists in stated directory
            self.assertTrue(os.path.exists(mod.solar_resource_file) ,
                            "Solar Resource file path could not be found for {0}".format(mod.__class__.__name__) )
    
if __name__ == "__main__":
    unittest.main()