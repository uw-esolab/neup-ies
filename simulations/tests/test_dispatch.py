#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:54:26 2021

@author: gabrielsoto
"""


from modules.NuclearTES import NuclearTES
from dispatch.NuclearDispatch import NuclearDispatch
import pyomo.environ as pe
from pyomo.environ import units as u_pyomo
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import unittest, os, math
import numpy as np

class TestGeneralDispatch(unittest.TestCase):
    """
    Unit tests for General Dispatch module
    
    This testing suite is meant to test the general Dispatch module of 
    the NE2 architecture. Here we have 
    """

    def setUp(self):
        """ Creating instance of GeneralDispatch upon start of each test
        """
        
        nuctes = NuclearTES(json_name='tests/test_nuctes',is_dispatch=True)
        
        # saving list of modules
        self.mod_list = [nuctes]
        
        # saving list of dispatch classes
        self.disp_list = [NuclearDispatch]
        
        # creating individual parameter dictionaries
        self.param_list = []
        
        for mod in self.mod_list:
            # ==============================================================
            # creating a duplicate Plant (with setup steps for time elements)
            
            # initialize times
            mod.run_loop = True
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
            
            # ==============================================================
            # create dispatch parameters dictionary
            params_dict = mod.create_dispatch_params( prePlant  )
            
            # log the dictionary
            self.param_list.append(params_dict)
            
            # delete extra new variablesmodel
            del time_start, time_next, prePlant, params_dict
            
    
    def tearDown(self):
        """ Deleting instance of GeneralDispatch at end of each test
        """

        # deleting each module
        for (mod, disp, params) in zip(self.mod_list, self.disp_list, self.param_list):
            del mod, disp, params
        
        # deleting module list. this might be overkill?
        del self.mod_list, self.disp_list, self.param_list


    def test__init__(self):
        """ Testing the shared processes in __init__ of all modules
        """

        param_attr_list = [ 'Ec',        # from GeneralDispatch.set_power_cycle_parameters()
                            'P',         # from GeneralDispatch.set_time_indexed_parameters()
                            'Cpc',       # from GeneralDispatch.set_fixed_cost_parameters()
                            'Cnuc',      # from NuclearDispatch.set_fixed_cost_parameters()
                            'En',        # from NuclearDispatch.set_nuclear_parameters()
                            'Qin_nuc',   # from NuclearDispatch.set_time_series_nuclear_parameters()
                            'y0',        # from GeneralDispatch.set_initial_state()
                            'yn0' ]      # from NuclearDispatch.set_initial_state()

        vars_attr_list = [  'wdot',      # from GeneralDispatch
                            'x',         # from GeneralDispatch
                            'xn',        # from NuclearDispatch
                            'yn']        # from NuclearDispatch
        
        # looping through all defined modules + attributes
        for (mod, disp, params) in zip(self.mod_list, self.disp_list, self.param_list):
            
            # ==============================================================
            # creating dispatch model
            disp_model = disp(params, mod.u)
            
            # ==============================================================
            # checking pyomo units and Pyomo Concrete model
            
            # check that pyomo units were created
            self.assertTrue(  hasattr(disp_model, 'u_pyomo'), "Pyomo units object not saved to dispatch class.")
            self.assertEqual( type(disp_model.u_pyomo ), type(u_pyomo), "u_pyomo not correct unit type.")
            
            # check that there exists a Pyomo model
            self.assertTrue( hasattr(disp_model, 'model'), "Pyomo model object not saved to dispatch class.")
            self.assertEqual( type(disp_model.model ), type(pe.ConcreteModel() ), "u_pyomo not correct unit type.")
            
            # ==============================================================
            # checking that parameters were saved to the Pyomo model
            
            # looping through all defined attributes
            for attr in param_attr_list:
                
                # checking that the attribute exists
                self.assertTrue( hasattr(disp_model.model, attr) ,
                            "Parameter {0} is not found in dispatch model for {1}".format(attr, disp) )

            # ==============================================================
            # checking that parameters were saved to the Pyomo model
            
            # looping through all defined attributes
            for attr in vars_attr_list:
                
                # checking that the attribute exists
                self.assertTrue( hasattr(disp_model.model, attr) ,
                            "Variable {0} is not found in dispatch model for {1}".format(attr, disp) )

            # ==============================================================
            # checking that objective function is saved to the Pyomo model
            
            self.assertTrue(  hasattr(disp_model.model, 'OBJ'), "Objective not saved to dispatch class.")
            

    def test_pyomo_model_units(self):
        """ Testing the unit consistency of dispatch models
        """
        
        # looping through all defined modules + attributes
        for (mod, disp, params) in zip(self.mod_list, self.disp_list, self.param_list):
            
            # ==============================================================
            # creating dispatch model
            disp_model = disp(params, mod.u)
            
            # ==============================================================
            # testing units of Objective
            assert_units_consistent( disp_model.model.OBJ )
    
            # ==============================================================
            # testing units of Constraints
            assert_units_consistent( disp_model.model )
            
        
if __name__ == "__main__":
    unittest.main()
            