#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:49:24 2021

@author: gabrielsoto
"""

import sys, copy
sys.path.append('..')
import numpy as np
import PySAM.TcsmoltenSalt as TcsmoltenSalt
from modules.GenericSSCModule import GenericSSCModule
from modules.NuclearTES import NuclearTES
from dispatch.SolarDispatch import SolarDispatch as SD
from dispatch.SolarDispatch import SolarDispatchParamWrap as SDP
from dispatch.SolarDispatch import SolarDispatchOutputs as SDO

class SolarTES(NuclearTES): 
    """
    The SolarTES class intializes, updates, and runs SSC simulations through PySAM,
    specifically for the SSC tcsmolten_salt module. 
    """
    
    def __init__(self, plant_name="tcsmolten_salt", json_name="model_solarTES", is_dispatch=False):
        """ Initializes the SolarTES module
        
        Inputs:
            plant_name (str)         : name of SSC module to run 
            json_name (str)          : name of JSON script with input data for module
            is_dispatch (bool)       : boolean, if True runs Pyomo dispatch optimization
        """
        
        # initialize Generic module, csv data arrays should be saved here
        NuclearTES.__init__( self, plant_name, json_name, is_dispatch )
        
        # define specific PySAM module to be called later
        self.PySAM_Module = TcsmoltenSalt
        
        # define specific Dispatch module to be called later
        self.Dispatch_Module = SD


        self.Dispatch_Outputs = SDO


    def update_Plant_after_SSC(self):
        """ Update SSC Plant inputs with SSC outputs from previous segment simulation
        
        ** self.run_loop == True
        
        This method uses the SSC end results from the previous simulation segment
        and sets them as the initial conditions for the next SSC segment. As a 
        small note: some outputs are arrays that span the full year, however the
        only relevant parts are the first indeces corresponding to the SSC Horizon.
        All other values are typically 0. 
        """
        
        ssc_slice = self.slice_ssc_firstH
        
        self.Plant.SystemControl.rec_op_mode_initial              = self.Plant.Outputs.rec_op_mode_final[ssc_slice][-1]
        self.Plant.SystemControl.rec_startup_time_remain_init     = self.Plant.Outputs.rec_startup_time_remain_final[ssc_slice][-1]
        self.Plant.SystemControl.rec_startup_energy_remain_init   = self.Plant.Outputs.rec_startup_energy_remain_final[ssc_slice][-1]
        self.Plant.SystemControl.T_tank_cold_init                 = self.Plant.Outputs.T_tes_cold[ssc_slice][-1]
        self.Plant.SystemControl.T_tank_hot_init                  = self.Plant.Outputs.T_tes_hot[ssc_slice][-1]
        self.Plant.ThermalStorage.csp_pt_tes_init_hot_htf_percent = self.Plant.Outputs.hot_tank_htf_percent_final[ssc_slice][-1]
        self.Plant.SystemControl.pc_op_mode_initial               = self.Plant.Outputs.pc_op_mode_final[ssc_slice][-1]
        self.Plant.SystemControl.pc_startup_energy_remain_initial = self.Plant.Outputs.pc_startup_time_remain_final[ssc_slice][-1]
        self.Plant.SystemControl.pc_startup_time_remain_init      = self.Plant.Outputs.pc_startup_energy_remain_final[ssc_slice][-1]
        self.Plant.SystemControl.is_field_tracking_init           = self.Plant.Outputs.is_field_tracking_final[ssc_slice][-1]


    def create_dispatch_wrapper(self, PySAM_dict):
        """ Creating a wrapper object for calling a class that creates dispatch parameters
        
        ** self.is_dispatch == True 
        (Called in __init__ of NE2 module)
        
        This method creates an object whose class ultimately calculates and creates 
        parameters for Dispatch optimization. The reason this class exists separately
        is that it gets overlaoded based on the PySAM module we are running. Depending on
        the PySAM module, this method calls on a different Dispatch Parameter class that 
        is specific to the module.

        Inputs:
            PySAM_dict (dict)   : dictionary of PySAM inputs from a script in /json directory
        Outputs:
            dispatch_wrap (obj) : wrapper object for the class that creates dispatch parameters
        """
        
        self.DispatchParameterClass = SDP
        
        dispatch_wrap = self.DispatchParameterClass( self.u, self.SSC_dict, PySAM_dict,
                    self.pyomo_horizon, self.dispatch_time_step)
        
        return dispatch_wrap


    def create_dispatch_params(self, Plant ):
        """ Populating a dictionary with dispatch parameters before optimization
        
        Note:
            self.is_dispatch == True 
            (Called within simulation)
        
        This method is creates the Dispatch Parameter dictionary that will be 
        populated with static inputs from SSC_dict as well as initial conditions
        for Dispatch optimization. The initial conditions are continuously updated
        if simulation is segmented.

        Args:
            Plant (obj): 
                original PySAM Plant module
        Returns:
            dispatch_wrap (obj): 
                wrapper object for the class that creates dispatch parameters
                
        """
        # get the object
        DW = self.dispatch_wrap
        
        # run the setters from the GenericSSCModule parent class
        params = GenericSSCModule.create_dispatch_params(self, Plant )
        
        # extract array from full run of Plant
        assert hasattr(Plant.Outputs, "Q_thermal"), "Q_thermal was not found in the outputs of Plant."
        self.Q_rec_guess = Plant.Outputs.Q_thermal
        
        # set up copy of SSC dict
        updated_SSC_dict = copy.deepcopy(self.SSC_dict)
        updated_SSC_dict['Q_thermal'] = self.Q_rec_guess[self.slice_pyo_firstH]
        
        # these are SolarTES-specific setters
        params = DW.set_solar_parameters( params )
        params = DW.set_time_series_solar_parameters( params, updated_SSC_dict )
        
        # this sets the initial states for the SolarTES
        params = DW.set_initial_state( params )
        
        return params
    
    
    def initialize_arrays(self):
        """ Initializing empty arrays to log SSC outputs after segment simulations
        
        This method creates empty arrays where SSC outputs will be written to.
        Also creates a list of str names for logged simulation outputs.
        
        """
        
        u = self.u
        
        # start and end times for full simulation
        i_start = (self.SSC_dict['time_start'] * u.s).to('hr').m
        i_end   = (self.SSC_dict['time_stop'] * u.s).to('hr').m
        
        # size of simulation arrays
        N_sim = int( i_end - i_start )
        
        # dictionary of output variable names to log after each segment simulation
        self.Log_Arrays = {
        #    name of NE2 variable || name of SSC module variable
                'time_log':          'time_hr',          # logging time
                'gen_log':           'gen',              # electricity generation log
                'q_thermal_log':     'Q_thermal',    # thermal power from nuclear to HTF 
                'p_cycle_log' :      'P_cycle',          # PC electrical power output (gross)
                'q_dot_rec_inc_log': 'q_dot_rec_inc',    # Nuclear incident thermal power
                'q_pb_log':          'q_pb',             # PC input energy
                'q_dot_pc_su_log' :  'q_dot_pc_startup', # PC startup thermal power
                'm_dot_pc_log' :     'm_dot_pc',         # PC HTF mass flow rate
                'm_dot_rec_log'  :   'm_dot_rec',        # Nuc mass flow rate
                'T_pc_in_log' :      'T_pc_in',          # PC HTF inlet temperature 
                'T_pc_out_log'   :   'T_pc_out',         # PC HTF outlet temperature
                'T_tes_cold_log':    'T_tes_cold',       # TES cold temperature
                'T_tes_hot_log'  :   'T_tes_hot',        # TES hot temperature
                'T_rec_in_log':      'T_rec_in',         # Plant inlet temperature
                'T_rec_out_log'  :   'T_rec_out',        # Plant outlet temperature
                'T_cond_out_log':    'T_cond_out',       # PC condenser water outlet temperature
                'e_ch_tes_log'  :    'e_ch_tes',         # TES charge state
                'op_mode_1_log' :    'op_mode_1',        # Operating Mode
                'defocus_log'   :    'defocus',          # Nuclear "Defocus" fraction
                'eta_log'       :    'eta'               # PC efficiency, gross
            } if self.run_loop \
                 else {'gen_log':    'gen'  # electricity generation log
                      }
        
        # empty array to initalize log arrays
        empty_array = np.zeros(N_sim)
        
        # loop through keys in ^ dictionary, save the KEY name to NE2 module as empty array
        for key in self.Log_Arrays.keys():
            # meta: if we don't grab the copy of empty_array, it'll assign a pointer to the array!!
            setattr( self, key, empty_array.copy() ) 
            
        if self.log_dispatch_targets:
            self.Log_Target_Arrays = {
                   'is_rec_su_allowed_in' : empty_array.copy(),
                   'is_rec_sb_allowed_in' : empty_array.copy(),
                   'is_pc_su_allowed_in'  : empty_array.copy(),
                   'is_pc_sb_allowed_in'  : empty_array.copy(),
                   'q_pc_target_su_in'    : empty_array.copy(),
                   'q_pc_target_on_in'    : empty_array.copy(),
                   'q_pc_max_in'          : empty_array.copy()
                   }

    def update_Pyomo_after_SSC(self, Plant, params ):
        """ Update Pyomo inputs with SSC outputs from previous segment simulation
        
        Note:
            self.run_loop    == True
            self.is_dispatch == True 
                          
        This method uses the SSC end results from the previous simulation segment
        and uses them to update the existing Dispatch parameter dictionary that
        is ultimately sent to Pyomo. Essentially just updates the initial conditions
        of the Dispatch parameter dictionary. 
        
        Args:
            Plant (obj): 
                original PySAM Plant module
            params (dict): 
                dictionary of Pyomo dispatch parameters
        Returns:
            params (dict): 
                updated dictionary of Pyomo dispatch parameters
                
        """
        
        ssc_slice = self.slice_ssc_firstH
        
        updated_SSC_dict = copy.deepcopy(self.SSC_dict)
        
        # saving relevant end-of-sim outputs from the last simulation segment
        updated_SSC_dict['rec_op_mode_initial']              = Plant.Outputs.rec_op_mode_final[ssc_slice][-1]
        updated_SSC_dict['rec_startup_time_remain_init']     = Plant.Outputs.rec_startup_time_remain_final[ssc_slice][-1]
        updated_SSC_dict['rec_startup_energy_remain_init']   = Plant.Outputs.rec_startup_energy_remain_final[ssc_slice][-1]
        updated_SSC_dict['T_tank_cold_init']                 = Plant.Outputs.T_tes_cold[ssc_slice][-1]
        updated_SSC_dict['T_tank_hot_init']                  = Plant.Outputs.T_tes_hot[ssc_slice][-1]
        updated_SSC_dict['csp.pt.tes.init_hot_htf_percent']  = Plant.Outputs.hot_tank_htf_percent_final[ssc_slice][-1]
        updated_SSC_dict['pc_op_mode_initial']               = Plant.Outputs.pc_op_mode_final[ssc_slice][-1]
        updated_SSC_dict['pc_startup_time_remain_init']      = Plant.Outputs.pc_startup_time_remain_final[ssc_slice][-1]
        updated_SSC_dict['pc_startup_energy_remain_initial'] = Plant.Outputs.pc_startup_energy_remain_final[ssc_slice][-1]
        updated_SSC_dict['is_field_tracking_init']           = Plant.Outputs.is_field_tracking_final[ssc_slice][-1]
        
        # these are specific to the initial states
        updated_SSC_dict['wdot0']     = Plant.Outputs.P_cycle[ssc_slice][-1]
        
        # extract time series from a previous SSC run
        updated_SSC_dict['Q_thermal'] = self.Q_rec_guess[self.slice_pyo_currentH]
        
        # TODO: removing w_dot_s_prev references in all of Dispatch for now, might need to revisit later
        # updated_SSC_dict['wdot_s_prev'] = 0 #np.array([pe.value(dm.model.wdot_s_prev[t]) for t in dm.model.T])[-1]
        
        DW = self.dispatch_wrap
        
        # set up Plant outputs dictionary
        plant_dict = Plant.Outputs.export()
        
        # updating the initial state and time series Nuclear params
        params = DW.set_time_indexed_parameters( params, self.df_array, self.ud_array, self.slice_pyo_currentH )
        params = DW.set_initial_state( params, updated_SSC_dict, Plant, self.t_ind )
        params = DW.set_time_series_solar_parameters( params, updated_SSC_dict )
        
        return params


    def update_Plant_after_Pyomo(self, Plant, pre_dispatch_run=False):
        """ Update SSC Plant inputs with Pyomo optimization outputs from current segment simulation

        Note:
            self.run_loop    == True (can be called outside loop)
            self.is_dispatch == True 
        
        This method uses the optimization results from Pyomo and ensures that 
        the next SSC segment uses them throughout the corresponding SSC Horizon.
        SSC normally takes single values for initial conditions (for the first hour
        of the SSC Horizon), but it can also take in an array of values for each
        step in the SSC Horizon. These are called "dispatch_targets". Steps are:
        (1) extract solutions from Pyomo over the Pyomo Horizon, 
        (2) keep the solutions for the shorter SSC Horizon and 
        (3) save these "dispatch target" inputs to the Plant object for the 
        next SSC simulation segment. 

        Args:
            Plant (obj): 
                original PySAM Plant module to be updated
            pre_dispatch_run (bool): 
                are we updating the Plant for a pre- or post- dispatch run.
                Recall that we only log post-dispatch Plant runs
        Returns:
            Plant (obj): 
                updated PySAM Plant module
                
        """
        
        # number of times in full simulation 
        N_full     = int((self.SSC_dict['time_stop']*self.u.s).to('hr').m)
        
        # the heavy-lifting happens here -> return a dictionary of dispatch target arrays from Pyomo optimization results
        horizon = self.pyomo_horizon if pre_dispatch_run else self.ssc_horizon
        dispatch_targets = self.Dispatch_Outputs.get_dispatch_targets_from_Pyomo( \
                                                self.current_disp_model, horizon, N_full, self.run_loop)
        
        ### Set Dispatch Targets ### 
        # setting dispatch targets to True so that SSC can read in Pyomo inputs
        Plant.SystemControl.is_dispatch_targets = True
        
        # extract binary arrays for receiver startup and standby
        Plant.SystemControl.is_rec_su_allowed_in = dispatch_targets['is_rec_su_allowed_in']
        Plant.SystemControl.is_rec_sb_allowed_in = dispatch_targets['is_rec_sb_allowed_in']
        
        # extract binary arrays for cycle startup and standby
        Plant.SystemControl.is_pc_su_allowed_in  = dispatch_targets['is_pc_su_allowed_in']
        Plant.SystemControl.is_pc_sb_allowed_in  = dispatch_targets['is_pc_sb_allowed_in']
        
        # extract power arrays for power cycle
        Plant.SystemControl.q_pc_target_su_in    = dispatch_targets['q_pc_target_su_in']
        Plant.SystemControl.q_pc_target_on_in    = dispatch_targets['q_pc_target_on_in']
        Plant.SystemControl.q_pc_max_in          = dispatch_targets['q_pc_max_in']
        
        return Plant