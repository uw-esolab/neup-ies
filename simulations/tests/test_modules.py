#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:36:02 2021

@author: gabrielsoto
"""

from modules.NuclearTES import NuclearTES
import unittest, os, math


class TestPySAMModules(unittest.TestCase):
    """
    Unit tests for PySAM module setup
    
    This testing suite is meant to test the methods in parent class of all modules,
    that being GenericSSCModule. Only non-abstract methods will be tested.
    All child classes will share common methods and attributes which will be tested here. 
    Dispatch is NOT tested here.
    
    For each individual class, there will be bespoke test classes.
    """

    def setUp(self):
        """ Creating instances of modules upon start of each test
        """
        
        nuctes = NuclearTES(json_name='tests/test_nuctes',is_dispatch=False)
        
        #saving list of modules
        self.mod_list = [nuctes]
    
    
    def tearDown(self):
        """ Deleting instances of modules at end of each test
        """
        
        # deleting each module
        for mod in self.mod_list:
            del mod
        
        # deleting module list. this might be overkill?
        del self.mod_list


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
            

    def test_initialize_arrays(self):
        """ Testing the initialize_arrays method
        """
        
        # running in time segments versus running full year at once
        loop_modes = [True, False]

        # looping through all modules
        for mod in self.mod_list:
            for mode in loop_modes:
                
                # =====================================================
                # run the method for given looping mode
                mod.run_loop = mode
                mod.initialize_arrays()
    
                # assert that a Log dictionary was created
                self.assertTrue( hasattr(mod, 'Log_Arrays') , 
                                "Log_Arrays dictionary not created in {0} module.".format(mod.__class__.__name__ ) )
    
                # check that individual log arrays were saved to the NE2 module
                keys = mod.Log_Arrays.keys()
                
                # if NOT running time segments, should only be saving the 'gen_long' array
                condition = len(keys) == 1 if mode is False else len(keys) > 1
                
                self.assertTrue( condition , 
                                "Improper amount of keys created in {0} module when run_loop is {1}.".format(mod.__class__.__name__ , mode ) )
                
                for k in keys:
                    self.assertTrue( hasattr(mod, k) , 
                                    "{0} Plant module doesn't have Output {1}.".format(mod.__class__.__name__ , k) )
                
                mod.reset_all()


    def test_reset_all(self):
        """ Testing the reset_all method
        """
        
        try_all_PySAM_mods = [True, False]
        
        # looping through all modules
        for mod in self.mod_list:
            for try_all in try_all_PySAM_mods:
                
                # create PySAM modules
                mod.create_Plant()
                
                # testing safe deletion, will it throw error if Grid was never created?
                if try_all:
                    mod.create_Grid()
                    mod.create_SO()
                
                # initialize arrays
                mod.run_loop = True
                mod.initialize_arrays()
                keys = mod.Log_Arrays.keys()
                
                # run the method
                mod.reset_all()
                
                # check that method deleted PySAM modules
                self.assertFalse( hasattr(mod,"Plant") ,
                                 "Plant not deleted in {0} module.".format(mod.__class__.__name__ ) )
                self.assertFalse( hasattr(mod,"Grid") ,
                                 "Grid not deleted in {0} module. Grid was {1} beforehand".format(mod.__class__.__name__ ,
                                                "created" if try_all else "NOT created") )
                self.assertFalse( hasattr(mod,"SO") ,
                                 "SO not deleted in {0} module. SO was {1} beforehand".format(mod.__class__.__name__ ,
                                                "created" if try_all else "NOT created") )
                
                # check that method deleted output and logging arrays
                self.assertFalse( hasattr(mod,"Log_Arrays") ,
                                 "Log_Arrays not deleted in {0} module.".format(mod.__class__.__name__ ) )
                
                for k in keys:
                    self.assertFalse( hasattr(mod, k) , 
                                    "{0} Plant module still has Output {1}".format(mod.__class__.__name__ , k) )
            
            
    # def test_log_SSC_arrays(self):
    #     """ Testing the log_SSC_arrays method
        
    #     Here we test both instances of the method call, during simulations and
    #     during the last simulation segment. 
    #     """
        
    #     # running in time segments versus running full year at once
    #     loop_modes = [True, False]
    #     final_attrs = ['gen_log', 'capacity_factor', 'annual_energy']
        
    #     # looping through all modules
    #     for mod in self.mod_list:
    #         # looping through all loop modes for log_final being true
    #         for mode in loop_modes: 
                
    #             # ================================================================
    #             # test that method behaves correctly for given loop mode
    #             mod.run_sim(run_loop=False)
                
                
    #             # check attributes 
    #             for f in final_attrs:
                    
    #                 # check that attributes exist
    #                 self.assertTrue( hasattr(mod, f) , 
    #                                 "{0} Plant module does not have output {1}".format(mod.__class__.__name__ , f) )
                    
    #                 # check that attribute is non-zero
    #                 self.assertTrue( getattr(mod, f).sum() != 0 ,
    #                                 "Output {0} sums up to 0 in module {1}".format(f, mod.__class__.__name__ ) )
                
    #             mod.reset_all()
            


    def test_run_sim(self):
        """ Testing run_sim for all modules
        
        NOTE: this doesn't test for accuracy of individual results, just that
        the processes run in the correct order and output something.
        """
        
        # list of important attributes for each submodule
        plant_output_attr = ['annual_energy', 'gen']
        grid_output_attr  = ['annual_energy_pre_curtailment_ac', 'gen']
        so_output_attr    = ['ppa', 'lppa_nom', 'lppa_real', 'project_return_aftertax_npv']
        
        # looping through all modules
        for mod in self.mod_list:
            
            # ======================================
            #---run full simulation for entire year
            mod.run_sim(run_loop=False)
            
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
            
            # ======================================
            #---run looped-simulation 
            # store results from full simulation
            annual_energy   = mod.Grid.SystemOutput.annual_energy
            ppa             = mod.SO.Outputs.ppa
            
            # reset submodules
            mod.reset_all()
            
            # ======================================
            #---run simulation in a loop
            mod.run_sim(run_loop=True)
            
            # check that results are in the same ballpark
            self.assertTrue( math.isclose (mod.Grid.SystemOutput.annual_energy, annual_energy, rel_tol=1e-2) , 
                             "Plant and Grid outputs are not within tolerance.")
            self.assertTrue( math.isclose (mod.SO.Outputs.ppa, ppa  , rel_tol=1e-2) , 
                             "Grid and SingleOwner outputs are not within tolerance.")
    
    
    def test_simulate_Plant(self):
        """ Testing the simulate_Plant method
        
        Here we test the basic functionality of simulate Plant
        """

        # looping through all modules
        for mod in self.mod_list:
            
            # =====================================================
            #---first test what happens when not running in loop
            mod.create_Plant( )
            mod.run_loop = False
            
            # simulate Plant
            mod.simulate_Plant( )
            
            # assert that sizing of time arrays are consistent
            t_ind     = mod.t_ind * mod.u.hr
            time_stop = mod.SSC_dict['time_stop'] * mod.u.s
            self.assertTrue( t_ind == time_stop.to('hr') , "Time index not set to time_stop when running full sim time." )
            
            # reset submodules
            del mod.Plant
            
            # =====================================================
            #---next test what happens when running in loop
            mod.create_Plant( )
            mod.run_loop = True
            
            # simulate Plant
            mod.simulate_Plant( )
            
            # assert that sizing of time arrays are consistent
            t_ind     = mod.t_ind * mod.u.hr
            time_stop = mod.ssc_horizon.to('s') 
            self.assertTrue( t_ind.to('s') == time_stop , "Time index not set to SSC Horizon time when running sim in loop." )


if __name__ == "__main__":
    unittest.main()