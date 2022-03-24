#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 14:39:17 2022

@author: gabrielsoto
"""


from modules.NuclearTES import NuclearTES
from dispatch.NuclearDispatch import NuclearDispatch
import unittest
import numpy as np


class TestNuclearTESScenarios(unittest.TestCase):
    """
    Unit tests for Nuclear TES module
    
    This testing suite is meant to test the nuclear TES module meant to
    represent Model 1 of the NE-2 project. This models an LFR plant supplying 
    power to a power block with thermal energy storage (TES) in the form 
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
    
    
    def test_dispatch_revenue(self):
        """ Testing that dispatch increases revenue and operations
        """
        
        # setting up modules
        mod_disp_true  = self.nuctes
        mod_disp_false = NuclearTES(json_name='tests/test_nuctes', is_dispatch=False)
        
        # running each module with and without pyomo dispatch optimization
        mod_disp_true.run_sim(run_loop=True)
        mod_disp_false.run_sim(run_loop=True)
        
        # determining total revenue from each module run
        revenue = {}
        mods = [mod_disp_true, mod_disp_false]
        names = ["w/ dispatch", 'w/o dispatch']
        for mod, name in zip(mods, names):
            
            # energy generated and associated price for each timestep
            gen   = np.array( mod.Plant.Outputs.gen )
            price = np.array( mod.Plant.TimeOfDeliveryFactors.dispatch_factors_ts )
            
            # logging revenue for each run
            revenue[name] = np.sum( gen * price )
        
        # asserting that revenue with dispatch is more than without dispatch
        self.assertTrue( revenue["w/ dispatch"] > revenue["w/o dispatch"], 
                          "Dispatch did not improve revenue for {0}".format( mod_disp_true ))
        
        
        dispatch_gen   = np.array( mod_disp_true.Plant.Outputs.gen )
        dispatch_price = np.array( mod_disp_true.Plant.TimeOfDeliveryFactors.dispatch_factors_ts )
        
        high_price = np.where( dispatch_price > 1 )[0]
        low_price  = np.where( dispatch_price < 1 )[0]
        
        # asserting that dispatch makes plant sell more energy during high price periods
        self.assertTrue( dispatch_gen[high_price].sum() > dispatch_gen[low_price].sum(), 
                          "Energy generation is not higher when prices are high for dispatch of {0}".format( mod_disp_true ) )


if __name__ == "__main__":
    unittest.main()