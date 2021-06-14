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
from dispatch.GeneralDispatch import GeneralDispatchOutputs as GDO
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
        
        # initialize dispatch-specific parameters and loggers
        if self.is_dispatch:
            
            # initialize dispatch counter
            self.disp_count = 0
            
            # empty dictionaries to store each dispatch run
            self.disp_models  = {}
            self.disp_results = {}
            
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
        
        if self.run_loop:
            self.log_SSC_arrays(log_final=True)
        
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
        self.t_ind = int(time_next.to('hr').m)
        
        # setting up empty log array for gen
        self.initialize_arrays()
        
        # run dispatch optimization for the first time
        if self.is_dispatch:
            # create dispatch parameters for the first time
            disp_params = self.create_dispatch_params()
            
            # run pyomo optimization
            self.run_pyomo(disp_params)
            
            # updating SSC inputs using Pyomo outputs
            self.update_Plant_after_Pyomo()
        
        # first execution of Plant through SSC
        self.run_Plant_through_SSC( time_start , time_next )
        
        # this loop should only be entered if run_loop == True
        while (time_next < time_end):
            
            # time-printer
            print_time = int(time_next.to('d').magnitude)
            if not print_time % 1: print('   [%s / %s] completed.' % (print_time, np.round(time_end.to('d').m)) )
            
            # update time
            time_start += self.ssc_horizon.to('s')
            time_next  += self.ssc_horizon.to('s')
            
            # update Plant parameters after previous run
            self.update_Plant_after_SSC( )
            
            # run dispatch optimization
            if self.is_dispatch:
                
                # i.e. updating Pyomo inputs using SSC outputs
                disp_params = self.update_Pyomo_after_SSC(disp_params)
                
                # run pyomo optimization again
                self.run_pyomo(disp_params)
            
                # updating SSC inputs using Pyomo outputs
                self.update_Plant_after_Pyomo( )
            
            # run Plant again
            self.run_Plant_through_SSC( time_start , time_next )
            
        
    def run_Plant_through_SSC(self, start_hr, end_hr):
        """ Simulation of Plant through SSC for given times
        """
        
        # start and end times for full simulation
        i_start = int( start_hr.to('hr').m )
        i_end   = int( end_hr.to('hr').m )
        
        self.Plant.SystemControl.time_start = start_hr.to('s').magnitude
        self.Plant.SystemControl.time_stop = end_hr.to('s').magnitude
        
        self.Plant.execute()
        
        # logging SSC outputs to arrays that have already been initialized
        if self.run_loop:
            self.log_SSC_arrays(i_start, i_end)
        else:
            self.gen_log[i_start:i_end] = self.Plant.Outputs.gen[0:self.t_ind ]
        
        
    def reset_all(self):
        """ Reset SSC submodules
        """
        del self.Plant
        del self.Grid
        del self.SO
        
    
    def run_pyomo(self, params):
        
        dispatch_model = GD(params, self.u)
        rt_results = dispatch_model.solve_model()
        
        self.current_disp_model = dispatch_model
        
        count = str(self.disp_count)
        self.disp_models[count]  = dispatch_model    
        self.disp_results[count] = rt_results
        
        self.disp_count += 1
    
    
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
      
        
    def update_Pyomo_after_SSC(self, params):
        return params
        
        
    def update_Plant_after_Pyomo(self):
        
        # number of times in full simulation 
        N_full     = int((self.SSC_dict['time_stop']*self.u.s).to('hr').m)
        
        # the heavy-lifting happens here -> return a dictionary of dispatch target arrays from Pyomo optimization results
        # NOTE: will need to overload GDO if this somehow gets overloaded by Nuclear Dispatch Outputs class
        dispatch_targets = GDO.get_dispatch_targets_from_Pyomo(self.current_disp_model, self.ssc_horizon, N_full, self.run_loop)
        
        ### Set Dispatch Targets ### 
        # setting dispatch targets to True so that SSC can read in Pyomo inputs
        self.Plant.SystemControl.is_dispatch_targets = True
        
        # extract binary arrays for receiver startup and standby
        self.Plant.SystemControl.is_rec_su_allowed_in = dispatch_targets['is_rec_su_allowed_in']
        self.Plant.SystemControl.is_rec_sb_allowed_in = dispatch_targets['is_rec_sb_allowed_in']
        
        # extract binary arrays for cycle startup and standby
        self.Plant.SystemControl.is_pc_su_allowed_in  = dispatch_targets['is_pc_su_allowed_in']
        self.Plant.SystemControl.is_pc_sb_allowed_in  = dispatch_targets['is_pc_sb_allowed_in']
        
        # extract power arrays for power cycle
        self.Plant.SystemControl.q_pc_target_su_in    = dispatch_targets['q_pc_target_su_in']
        self.Plant.SystemControl.q_pc_target_on_in    = dispatch_targets['q_pc_target_on_in']
        self.Plant.SystemControl.q_pc_max_in          = dispatch_targets['q_pc_max_in']
        

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

        return params

    def initialize_arrays(self):
        
        u = self.u
        
        # start and end times for full simulation
        i_start = (self.SSC_dict['time_start'] * u.s).to('hr').m
        i_end   = (self.SSC_dict['time_stop'] * u.s).to('hr').m
        
        # size of simulation arrays
        N_sim = int( i_end - i_start )
        
        # setting each logging array to zero
        self.time_log          = np.zeros(N_sim)
        self.gen_log           = np.zeros(N_sim)
        self.p_cycle_log       = np.zeros(N_sim)
        self.q_dot_rec_inc_log = np.zeros(N_sim)
        self.m_dot_pc_log      = np.zeros(N_sim)
        self.m_dot_rec_log     = np.zeros(N_sim)
        self.T_pc_in_log       = np.zeros(N_sim)
        self.T_pc_out_log      = np.zeros(N_sim)
        self.e_ch_tes_log      = np.zeros(N_sim)
        self.op_mode_1_log     = np.zeros(N_sim)
        self.defocus_log       = np.zeros(N_sim)


    def log_SSC_arrays(self, i_start=0, i_end=1, log_final=False):
        
        ssch   = slice(i_start,i_end,1)
        firsth = slice(0,self.t_ind,1)
        
        if not log_final:
            self.time_log[ssch]          = self.Plant.Outputs.time_hr[firsth]
            self.gen_log[ssch]           = self.Plant.Outputs.gen[firsth]
            self.p_cycle_log[ssch]       = self.Plant.Outputs.P_cycle[firsth]
            self.q_dot_rec_inc_log[ssch] = self.Plant.Outputs.q_dot_rec_inc[firsth]
            self.m_dot_pc_log[ssch]      = self.Plant.Outputs.m_dot_pc[firsth]
            self.m_dot_rec_log[ssch]     = self.Plant.Outputs.m_dot_rec[firsth]
            self.T_pc_in_log[ssch]       = self.Plant.Outputs.T_pc_in[firsth]
            self.T_pc_out_log[ssch]      = self.Plant.Outputs.T_pc_out[firsth]
            self.e_ch_tes_log[ssch]      = self.Plant.Outputs.e_ch_tes[firsth]
            self.op_mode_1_log[ssch]     = self.Plant.Outputs.op_mode_1[firsth]
            self.defocus_log[ssch]       = self.Plant.Outputs.defocus[firsth]

        else:
            # wanted to create a quick subclass that where I can extract things during PostProcessing steps...
            self.Plant.PySAM_Outputs = lambda: None # don't try this at home
            
            convert_output = lambda arr: tuple( arr.tolist() )
            
            self.Plant.PySAM_Outputs.time_hr       = convert_output( self.time_log )
            self.Plant.PySAM_Outputs.gen           = convert_output( self.gen_log )
            self.Plant.PySAM_Outputs.P_cycle       = convert_output( self.p_cycle_log )
            self.Plant.PySAM_Outputs.q_dot_rec_inc = convert_output( self.q_dot_rec_inc_log )
            self.Plant.PySAM_Outputs.m_dot_pc      = convert_output( self.m_dot_pc_log )
            self.Plant.PySAM_Outputs.m_dot_rec     = convert_output( self.m_dot_rec_log )
            self.Plant.PySAM_Outputs.T_pc_in       = convert_output( self.T_pc_in_log )
            self.Plant.PySAM_Outputs.T_pc_out      = convert_output( self.T_pc_out_log )
            self.Plant.PySAM_Outputs.e_ch_tes      = convert_output( self.e_ch_tes_log )
            self.Plant.PySAM_Outputs.op_mode_1     = convert_output( self.op_mode_1_log )
            self.Plant.PySAM_Outputs.defocus       = convert_output( self.defocus_log )

    
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

            