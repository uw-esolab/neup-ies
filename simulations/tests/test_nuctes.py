#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

from modules.NuclearTES import NuclearTES
from dispatch.NuclearDispatch import NuclearDispatch
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
        self.nuctes = NuclearTES(json_name='tests/test_nuctes', is_dispatch=True)
        self.nuctes_name = self.nuctes.__class__.__name__
        
    
    def tearDown(self):
        """ Deleting instance oself.nuctesf NuclearTES at end of each test
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
                 
        
        mod = self.nuctes
        
        # ==============================================================
        # creating a duplicate Plant (with setup steps for time elements)
        
        # initialize times
        mod.run_loop = False
        time_start, time_next = mod.initialize_time_elements()
        mod.initialize_time_slices( time_start )
        mod.initialize_arrays()
        
        # create Plant
        mod.create_Plant()
        prePlant = mod.duplicate_Plant( mod.Plant )

        # ==============================================================
        # pre-run the duplicate plant to gather guesses for Q_nuc
        
        ssc_run_success, prePlant = mod.run_Plant_through_SSC(
                                            prePlant, time_start , mod.sim_time_end 
                                            )
        self.assertTrue( ssc_run_success , 
                            "Pre-run of simulation failed for {0}".format( mod ) )
        
        # ==============================================================
        # actually run the create_dispatch_params method
        
        # create the dictionary
        params_dict = mod.create_dispatch_params( prePlant )
        
        # check that pre-run worked
        self.assertTrue( hasattr(mod, 'Q_nuc_guess'), 
                          "Guess for nuclear q_dot profile is not saved to {0} from simulation pre-run".format( mod ))

        self.assertTrue( sum( mod.Q_nuc_guess ) > 0, 
                          "Nuclear q_dot profile values are all 0 for {0}".format( mod ))
        
        # looping through all defined attributes
        for attr in attr_list:
            
            # checking that the attribute exists
            self.assertIn( attr, params_dict.keys() ,
                            "Parameter {0} is not found in dispatch parameter dictionary for {1}".format(attr, mod) )
            
        # erase mod
        del mod
        

if __name__ == "__main__":
    unittest.main()
            