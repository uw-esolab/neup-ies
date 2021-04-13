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
        genmod = GenericSSCModule()
        nuctes = NuclearTES()
        
        #saving list of modules
        self.mod_list = [genmod, nuctes]
        self.mod_list_CL = [nuctes]
    
    
    def tearDown(self):
        """ Deleting instances of modules at end of each test
        """
        
        # deleting each module
        for mod in self.mod_list:
            del mod
        
        for mod in self.mod_list_CL:
            del mod
        
        # deleting module list. this might be overkill?
        del self.mod_list
        del self.mod_list_CL


    def test__init__(self):
        """ Testing the shared processes in __init__ of all modules
        """
        
        # defining necessary attributes to be declared in __init__
        attr_list = ['json_name' , 'plant_name', 'ssc_horizon', 'pyomo_horizon', 
                     'SSC_dict']
        
        # looping through all defined modules + attributes
        for mod in self.mod_list:
            for attr in attr_list:
                
                # checking that all attributes exist in every module
                self.assertTrue(hasattr(mod,attr) , 
                                "{0} does not have '{1}' attribute".format(mod.__class__.__name__,attr))
            
            # checking that the self.store_csv_arrays( ) method was called correctly
            self.assertTrue(hasattr(mod,"solar_resource_file") ,
                            "Something went wrong when {0} called 'store_csv_arrays' method".format(mod.__class__.__name__) )


    def test_store_csv_arrays(self):
        """ Testing the storage of csv arrays in all modules
        
        NOTE: this assumes that the method was already called in __init__
        """
        
        # looping through all modules
        for mod in self.mod_list:
            
            # checking that the solar_resource_file path attribute exists
            self.assertTrue(hasattr(mod,"solar_resource_file") ,
                            "Something went wrong when {0} called 'store_csv_arrays' method".format(mod.__class__.__name__) )
            
            # checking that the actual filepath exists in stated directory
            self.assertTrue(os.path.exists(mod.solar_resource_file) ,
                            "Solar Resource file path could not be found for {0}".format(mod.__class__.__name__) )
    

    def test_create_Plant(self):
        """ Testing the create_Plant method
        """
        
        # looping through all modules
        for mod in self.mod_list:
            
            mod.create_Plant()
            # checking that Plant object was created
            self.assertTrue(hasattr(mod,'Plant') ,
                            "Plant object not created for {0}".format(mod.__class__.__name__) )
   
            
    def test_create_Grid(self):
        """ Testing the create_Grid method
        """
        
        # looping through all modules
        for mod in self.mod_list:
            
            mod.create_Plant()
            mod.create_Grid()
            
            # checking that Grid object was created
            self.assertTrue(hasattr(mod,'Grid') ,
                            "Grid object not created for {0}".format(mod.__class__.__name__) )


    def test_create_SO(self):
        """ Testing the create_SO method
        """
        
        # looping through all modules
        for mod in self.mod_list:
            
            mod.create_Plant()
            mod.create_SO()
            
            # checking that SingleOwner object was created
            self.assertTrue(hasattr(mod,'SO') ,
                            "SO object not created for {0}".format(mod.__class__.__name__) )
            
    
    def test_run_sim(self):
        """ Testing run_sim for all modules
        
        NOTE: this doesn't test for accuracy of individual results, just that
        the processes run in the correct order and output something.
        """
        
        plant_output_attr = ['annual_energy', 'gen']
        grid_output_attr  = ['annual_energy_pre_curtailment_ac', 'gen']
        so_output_attr    = ['ppa', 'lppa_nom', 'lppa_real', 'project_return_aftertax_npv']
        
        # looping through all modules
        for mod in self.mod_list:
            
            # run full simulation for entire year
            mod.run_sim(0)
            
            # check that Plant have written outputs
            for p_attr in plant_output_attr:
                self.assertTrue(hasattr(mod.Plant.Outputs, p_attr) ,
                                "{0} Plant doesn't have Output {1}".format(mod.__class__.__name__ , p_attr) )
            
            # check that Grid have written outputs
            for g_attr in grid_output_attr:
                self.assertTrue(hasattr(mod.Grid.Outputs, g_attr) ,
                                "{0} Grid doesn't have Output {1}".format(mod.__class__.__name__ , g_attr) )
            
            # check that SO have written outputs
            for s_attr in so_output_attr:
                self.assertTrue(hasattr(mod.SO.Outputs, s_attr) ,
                                "{0} SO doesn't have Output {1}".format(mod.__class__.__name__ , s_attr) )
            
            # fs = copy.deepcopy(mod)
            
            # if mod in self.mod_list_CL:
            #     mod.run_sim(1)
        
    
        
if __name__ == "__main__":
    unittest.main()