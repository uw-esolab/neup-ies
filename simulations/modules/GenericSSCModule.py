#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:40:22 2021

@author: gabrielsoto
"""

import PySAM.GenericSystem as GenericSystem
import PySAM.Grid as Grid
import PySAM.Singleowner as Singleowner
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
import pint
u = pint.UnitRegistry()
import numpy as np
import copy, os

class GenericSSCModule(object):
    
    def __init__(self, plant_name="generic_system", json_name="100mW_Generic"):
        
        # grab names, either default here or from child class
        self.json_name  = json_name
        self.plant_name = plant_name
        
        # read in dictionaries from json script
        PySAM_dict, SSC_dict = FileMethods.read_json( self.json_name )
        
        # storing SSC and Pyomo time horizons, inputs are in unit of hours
        self.ssc_horizon   = PySAM_dict['ssc_horizon'] * u.hr
        self.pyomo_horizon = PySAM_dict['pyomo_horizon'] * u.hr
        
        # save SSC_dict for usage later
        self.SSC_dict = SSC_dict
        
        # save csv arrays to class 
        self.store_csv_arrays( PySAM_dict )


    def run_sim(self, run_loop=False):
        """ Method to run single simulation for Generic System
        """
        
        self.run_loop = run_loop
        SSC_dict = self.SSC_dict
        
        #--- create Plant object and execute it
        self.create_Plant( )
        self.simulate_Plant( )
        
        #logging annual energy and capacity factor
        annual_energy = np.sum(self.gen_log)*u.kWh
        if 'P_ref' in SSC_dict.keys():
            ref_energy    = SSC_dict['P_ref']*u.MW * SSC_dict['gross_net_conversion_factor']
            self.capacity_factor = (annual_energy / (ref_energy * 1*u.yr) ).to(' ')
        
        #--- use executed Plant object to create Grid object and execute it
        self.create_Grid( )
        self.Grid.SystemOutput.gen = tuple(self.gen_log)
        self.Grid.SystemOutput.annual_energy = np.sum(annual_energy.magnitude)
        self.Grid.execute( )
        
        #--- use executed Plant object to create SingleOwner object and execute it
        self.create_SO( )
        self.SO.execute( )
        
         
    def store_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files and store in class
        
        Inputs:
            input_dict (dict) : dictionary with csv file names
        """
        
        # saving location of solar resource file for SSC input
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = os.path.join(parent_dir, input_dict['solar_resource_rel_parent']) #os.path.join
        
        
    def create_Plant(self):
        """ Method to create Plant object for the first time
        """
        
        # create plant data encoding for generic system
        plant_dat = pssc.dict_to_ssc_table( self.SSC_dict, self.plant_name )
        
        # create new Plant object
        self.Plant = GenericSystem.wrap(plant_dat)

    
    def create_Grid(self):
        """ Method to create Grid object for the first time
        """
        
        # create grid data encoding for grid
        grid_dat = pssc.dict_to_ssc_table( self.SSC_dict, "grid" )
        
        # create new Grid object from existing Plant object
        self.Grid = Grid.from_existing( self.Plant )
        
        # import Grid-specific data to Grid object
        self.Grid.assign(Grid.wrap(grid_dat).export())


    def create_SO(self):
        """ Method to create SingleOwner object for the first time
        """
        
        # create singleowner data encoding for singleowner object
        so_dat   = pssc.dict_to_ssc_table( self.SSC_dict, "singleowner" )
        
        # create new Singleowner object from existing Plant object
        self.SO = Singleowner.from_existing( self.Plant )
        
        # import Singleowner-specific data to Singleowner object
        self.SO.assign(Singleowner.wrap(so_dat).export())


    def simulate_Plant(self):
        """ Full simulation of Plant
        """
        
        # helper method for hstack-ing arrays
        augment_log = lambda X,Y: np.hstack(  [ X, Y ]  )
        
        # start and end times for full simulation
        time_start = self.SSC_dict['time_start'] * u.s
        time_end   = self.SSC_dict['time_stop'] * u.s
        
        # if running loop -> next time stop is ssc_horizon time
        #            else -> next time stop is sim end time
        time_next  = self.ssc_horizon.to('s') if self.run_loop else copy.deepcopy(time_end)
        
        # setting index for subsequent calls, static index for gen
        self.t_ind = int(time_next.to('hr').magnitude)
        
        # first execution of Plant through SSC
        self.run_Plant_through_SSC( time_start , time_next )
        
        # logging values of gen
        self.gen_log = np.ndarray([0])
        self.gen_log = augment_log( self.gen_log, self.Plant.Outputs.gen[0:self.t_ind ] )
        
        # this loop should only be entered if run_loop == True
        while (time_next < time_end):
            print_time = int(time_next.to('d').magnitude)
            if not print_time % 50: print('   [%s / %s] completed.' % (print_time,time_end.to('d').magnitude))
            
            # update time
            time_start += self.ssc_horizon.to('s')
            time_next  += self.ssc_horizon.to('s')
            
            # update Plant parameters after previous run
            self.update_Plant( )
            
            # run Plant again
            self.run_Plant_through_SSC( time_start , time_next )
            
            # log results
            self.gen_log = augment_log( self.gen_log, self.Plant.Outputs.gen[0:self.t_ind ] )
            

    def run_Plant_through_SSC(self, start_hr, end_hr):
        """ Simulation of Plant through SSC for given times
        """
        
        # oops, GenericSystem doesn't have 'SystemControl' subclass...
        if hasattr(self.Plant,'SystemControl'):
            self.Plant.SystemControl.time_start = start_hr.to('s').magnitude
            self.Plant.SystemControl.time_stop = end_hr.to('s').magnitude
        
        self.Plant.execute()
        
        
    def reset_all(self):
        """ Reset SSC submodules
        """
        del self.Plant
        del self.Grid
        del self.SO
        
    
    
    def update_Plant(self):
        """ Replace key Plant inputs with previous outputs in chainlinked simulations
        """
        
        self.Plant.SystemControl.rec_op_mode_initial              = self.Plant.Outputs.rec_op_mode_final
        self.Plant.SystemControl.rec_startup_time_remain_init     = self.Plant.Outputs.rec_startup_time_remain_final
        self.Plant.SystemControl.rec_startup_energy_remain_init   = self.Plant.Outputs.rec_startup_energy_remain_final
        self.Plant.SystemControl.T_tank_cold_init                 = self.Plant.Outputs.T_tes_cold[self.t_ind-1]
        self.Plant.SystemControl.T_tank_hot_init                  = self.Plant.Outputs.T_tes_hot[self.t_ind-1]
        self.Plant.ThermalStorage.csp_pt_tes_init_hot_htf_percent = self.Plant.Outputs.hot_tank_htf_percent_final
        self.Plant.SystemControl.pc_op_mode_initial               = self.Plant.Outputs.pc_op_mode_final
        self.Plant.SystemControl.pc_startup_energy_remain_initial = self.Plant.Outputs.pc_startup_time_remain_final
        self.Plant.SystemControl.pc_startup_time_remain_init      = self.Plant.Outputs.pc_startup_energy_remain_final
        
        