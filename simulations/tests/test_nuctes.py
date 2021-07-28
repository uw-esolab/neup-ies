#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

from modules.NuclearTES import NuclearTES
from dispatch.NuclearDispatch import NuclearDispatch
import unittest
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent


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
        self.nuctes = NuclearTES(json_name='tests/test_nuctes', is_dispatch=True)
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


    def test_create_dispatch_params(self):
        """ Testing the creation of a parameter dictionary for NuclearDispatch
        """
        
        attr_list = ['Ec',        # from GeneralDispatch.set_power_cycle_parameters()
                     'P',         # from GeneralDispatch.set_time_indexed_parameters()
                     'Cpc',       # from GeneralDispatch.set_fixed_cost_parameters()
                     'Cnuc',      # from NuclearDispatch.set_fixed_cost_parameters()
                     'En',        # from NuclearDispatch.set_nuclear_parameters()
                     'Qin_nuc',   # from NuclearDispatch.set_time_series_nuclear_parameters()
                     'y0',        # from GeneralDispatch.set_initial_state()
                     'yn0' ]      # from NuclearDispatch.set_initial_state()
                 
        
        # feigning a pyomo horizon slice
        time_slice = slice(0,48,1)
        
        # create the dictionary
        params_dict = self.nuctes.create_dispatch_params( time_slice )
        
        # looping through all defined attributes
        for attr in attr_list:
            
            # checking that the attribute exists
            self.assertIn( attr, params_dict.keys() ,
                            "Parameter {0} is not found in dispatch parameter dictionary for {1}".format(attr, self.nuctes_name) )


    def test_pyomo_model_units(self):
        """ Testing the unit consistency of NuclearDispatch model
        """
        
        # feigning a pyomo horizon slice
        time_slice = slice(0,48,1)
        
        # create the dictionary
        params_dict = self.nuctes.create_dispatch_params( time_slice )
        
        dispatch_model = NuclearDispatch( params_dict, self.nuctes.u )
        
        # actual pyomo model
        model = dispatch_model.model
        u_pyomo = dispatch_model.u_pyomo
        
        # ==============================================================
        # testing units of Objective
        
        assert_units_consistent( dispatch_model.model.OBJ )

        # ==============================================================
        # testing units of Constraints
        
        assert_units_consistent( dispatch_model.model )
        

if __name__ == "__main__":
    unittest.main()
            