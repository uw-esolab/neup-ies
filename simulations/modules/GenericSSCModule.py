#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:40:22 2021

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import PySAM.GenericSystem as GenericSystem
import PySAM.Grid as Grid
import PySAM.Singleowner as Singleowner
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
from dispatch.GeneralDispatch import GeneralDispatch as GD
from dispatch.GeneralDispatch import GeneralDispatchParamWrap as GDP
import pyomo.environ as pe
import numpy as np
import copy
from abc import ABC

class GenericSSCModule(ABC):
    
    def __init__(self, plant_name="nuclear_tes", json_name="model1", 
                       is_dispatch=False, dispatch_time_step=1):
        
        # grab names, either default here or from child class
        self.json_name  = json_name
        self.plant_name = plant_name
        
        # create and save a specific unit registry
        self.u = SSCHelperMethods.define_unit_registry()
        u = self.u
        
        # read in dictionaries from json script
        PySAM_dict, SSC_dict, output_keys = FileMethods.read_json( self.json_name )
        
        # storing SSC and Pyomo time horizons, inputs are in unit of hours
        self.ssc_horizon   = PySAM_dict['ssc_horizon'] * u.hr
        self.pyomo_horizon = PySAM_dict['pyomo_horizon'] * u.hr
        self.output_keys   = output_keys
        self.dispatch_time_step = dispatch_time_step * u.hr
        
        # save SSC_dict for usage later
        self.SSC_dict = SSC_dict
        
        # save csv arrays to class 
        self.store_csv_arrays( PySAM_dict )
        
        # save flag for dispatch 
        self.is_dispatch = is_dispatch
        
        if self.is_dispatch:
            # initialize dispatch wrap class
            self.dispatch_wrap = self.create_dispatch_wrapper( PySAM_dict )


    def run_sim(self, run_loop=False, export=False, filename='temp.csv'):
        """ Method to run single simulation for Generic System
        """
        
        u = self.u
        
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
        # update gen and annual energy so SystemOutput is consistent, carried over to SO object
        self.Grid.SystemOutput.gen = tuple(self.gen_log)
        self.Grid.SystemOutput.annual_energy = np.sum(annual_energy.magnitude)
        self.Grid.execute( )
        
        #--- use executed Plant object to create SingleOwner object and execute it
        self.create_SO( )
        self.SO.execute( )
        
        if export:
            self.export_results(filename)
         
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
        self.augment_log = lambda X,Y: np.hstack(  [ X, Y ]  )
        
        u = self.u
        
        # start and end times for full simulation
        time_start = self.SSC_dict['time_start'] * u.s
        time_end   = self.SSC_dict['time_stop'] * u.s
        
        # if running loop -> next 'time stop' is ssc_horizon time
        #            else -> next 'time stop' is sim end time
        time_next  = self.ssc_horizon.to('s') if self.run_loop else copy.deepcopy(time_end)
        
        # setting index for subsequent calls, static index for gen
        self.t_ind = int(time_next.to('hr').magnitude)
        
        # setting up empty log array for gen
        self.gen_log = np.ndarray([0])
        
        # run dispatch optimization for the first time
        if self.is_dispatch:
            
            # create dispatch parameters for the first time
            disp_params = self.create_dispatch_params()
            
            # run pyomo optimization
            outputs = self.run_pyomo(disp_params)
            
            # updating SSC inputs using Pyomo outputs
            self.update_Plant_after_Pyomo()
        
        # first execution of Plant through SSC
        self.run_Plant_through_SSC( time_start , time_next )
        
        
        # this loop should only be entered if run_loop == True
        while (time_next < time_end):
            print_time = int(time_next.to('d').magnitude)
            if not print_time % 50: print('   [%s / %s] completed.' % (print_time,time_end.to('d').magnitude))
            
            # update time
            time_start += self.ssc_horizon.to('s')
            time_next  += self.ssc_horizon.to('s')
            
            # update Plant parameters after previous run
            self.update_Plant_after_SSC( )
            
            # run dispatch optimization
            if self.is_dispatch:
                pass
                # **update** dispatch parameters using previous SSC run
                # i.e. updating Pyomo inputs using SSC outputs
                # disp_params = self.create_dispatch_params()
                
                # run pyomo optimization again
                # outputs = self.run_pyomo(disp_params)
            
                # # updating SSC inputs using Pyomo outputs
                # self.update_Plant_after_Pyomo( )
            

            # run Plant again
            self.run_Plant_through_SSC( time_start , time_next )
            
        
    def run_Plant_through_SSC(self, start_hr, end_hr):
        """ Simulation of Plant through SSC for given times
        """
        
        self.Plant.SystemControl.time_start = start_hr.to('s').magnitude
        self.Plant.SystemControl.time_stop = end_hr.to('s').magnitude
        
        self.Plant.execute()
        
        # logging values of gen
        self.gen_log = self.augment_log( self.gen_log, self.Plant.Outputs.gen[0:self.t_ind ] )
        
        
    def reset_all(self):
        """ Reset SSC submodules
        """
        del self.Plant
        del self.Grid
        del self.SO
        
    
    def run_pyomo(self, params):
        
        dispatch_model = GD(params, self.u)
        rt_results = dispatch_model.solve_model()
        
        self.dispatch_model = dispatch_model
        self.rt_results = rt_results
        return rt_results
    
    
    def update_Plant_after_SSC(self):
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
        
        
    def update_Plant_after_Pyomo(self):
        
        
        # setting dispatch targets to True so that SSC can read in Pyomo inputs
        self.Plant.SystemControl.is_dispatch_targets = True
        
        dm = self.dispatch_model
        
        # TODO: I think all of this stuff should live in GeneralDispatch or the Wrapper
        
        t_range = dm.model.T
        N_dispatch = dm.model.num_periods.value
        N_full     = int((self.SSC_dict['time_stop']*self.u.s).to('hr').m)
        empty_array = [0]*N_full
        horizon = slice(0,self.pyomo_horizon.m,1)
        # TODO: ^ make full 8760 sized array with front 48 spaces filled and save it to Plant
        
        # rec stuff
        yr   = np.array([pe.value(dm.model.yr[t]) for t in t_range])
        yrsu = np.array([pe.value(dm.model.yrsu[t]) for t in t_range])
        yrsb = np.array([pe.value(dm.model.yrsb[t]) for t in t_range])
        
        # cycle stuff
        x    = np.array([pe.value(dm.model.x[t]) for t in t_range])/1000. # from kWt -> MWt
        y    = np.array([pe.value(dm.model.y[t]) for t in t_range])
        ycsu = np.array([pe.value(dm.model.ycsu[t]) for t in t_range])
        ycsb = np.array([pe.value(dm.model.ycsb[t]) for t in t_range])
        
        Qc = np.array([pe.value(dm.model.Qc[t]) for t in t_range])/1000. # from kWt -> MWt
        Qu = dm.model.Qu.value/1000. # from kWt -> MWt

        # save partial arrays
        is_rec_su_allowed_in = [ 1 if (yr[t] + yrsu[t] + yrsb[t]) > 0.001 else 0 for t in range(N_dispatch)]  # Receiver on, startup, or standby
        is_rec_sb_allowed_in = [ 1 if yrsb[t] > 0.001 else 0 for t in range(N_dispatch)]  # Receiver standby

        is_pc_su_allowed_in = [ 1 if (y[t] + ycsu[t]) > 0.001 else 0 for t in range(N_dispatch)]  # Cycle on or startup
        is_pc_sb_allowed_in = [ 1 if ycsb[t] > 0.001 else 0 for t in range(N_dispatch)]  # Cyle standby

        #TODO: Might need to modify q_pc_target_on_in and q_pc_max_in for timesteps split between cycle startup and operation (e.g. 1383 - 1414 of csp_solver_core.cpp in mjwagner2/ssc/daotk-develop)
        q_pc_target_su_in = [Qc[t] if ycsu[t] > 0.001 else 0.0 for t in range(N_dispatch)]
        q_pc_target_on_in = [x[t] for t in range(N_dispatch)]
        q_pc_max_in = [Qu for t in range(N_dispatch)]
        
        
        # save full arrays
        self.Plant.SystemControl.is_rec_su_allowed_in = [ 1 if (yr[t] + yrsu[t] + yrsb[t]) > 0.001 else 0 for t in range(N_dispatch)]  # Receiver on, startup, or standby
        self.Plant.SystemControl.is_rec_sb_allowed_in = [ 1 if yrsb[t] > 0.001 else 0 for t in range(N_dispatch)]  # Receiver standby

        self.Plant.SystemControl.is_pc_su_allowed_in = [ 1 if (y[t] + ycsu[t]) > 0.001 else 0 for t in range(N_dispatch)]  # Cycle on or startup
        self.Plant.SystemControl.is_pc_sb_allowed_in = [ 1 if ycsb[t] > 0.001 else 0 for t in range(N_dispatch)]  # Cyle standby

        #TODO: Might need to modify q_pc_target_on_in and q_pc_max_in for timesteps split between cycle startup and operation (e.g. 1383 - 1414 of csp_solver_core.cpp in mjwagner2/ssc/daotk-develop)
        self.Plant.SystemControl.q_pc_target_su_in = [Qc[t] if ycsu[t] > 0.001 else 0.0 for t in range(N_dispatch)]
        self.Plant.SystemControl.q_pc_target_on_in = [x[t] for t in range(N_dispatch)]
        self.Plant.SystemControl.q_pc_max_in = [Qu for t in range(N_dispatch)]


    def create_dispatch_wrapper(self, PySAM_dict):
        
        self.DispatchParameterClass = GDP
        
        dispatch_wrap = self.DispatchParameterClass( self.u, self.SSC_dict, PySAM_dict,
                    self.pyomo_horizon, self.dispatch_time_step)
        
        return dispatch_wrap


    def create_dispatch_params(self):
        
        params = {}
        DW = self.dispatch_wrap
        
        # setting parameters for the first time
        params = DW.set_time_indexed_parameters( params )
        params = DW.set_power_cycle_parameters( params, self.ud_array )
        params = DW.set_fixed_cost_parameters( params )

        ### Initial Condition Parameters ###
        
        
        return params
        
    def export_results(self, filename):
        
        location_list = []
        key_list      = []
        value_list    = []
        
        for key in self.output_keys['keywords']:
            if hasattr(self,key):
                location_list.append('PySAM')
                value_list.append(getattr(self,key))
                
            elif hasattr(self.Plant.Outputs,key):
                location_list.append('Plant')
                value_list.append(getattr(self.Plant.Outputs,key))
                
            elif hasattr(self.SO.Outputs,key):
                location_list.append('SingleOwner')
                value_list.append(getattr(self.SO.Outputs,key))
                
            else:
                location_list.append('Not Found')
                value_list.append(' ')
                
            key_list.append(key)


        dataframe_list = [[x,y,z] for x,y,z in zip(location_list,key_list,value_list)]
        columns=['Object', 'Key', 'Value']
        
        FileMethods.write_csv(dataframe_list, columns, filename)

            