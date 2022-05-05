#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 16:58:30 2022

@author: gabrielsoto
"""

import sys, copy, os
sys.path.append('..')
import numpy as np
from util.FileMethods import FileMethods
from scipy.interpolate import interp1d
from util.SSCHelperMethods import SSCHelperMethods
import PySAM.NuclearMsptIndirectTes as NuclearMsptIndirectTes
from modules.GenericSSCModule import GenericSSCModule
from modules.NuclearTES import NuclearTES
from dispatch.IndirectNuclearDispatch import IndirectNuclearDispatch as IND
from dispatch.IndirectNuclearDispatch import IndirectNuclearDispatchParamWrap as INDP
from dispatch.IndirectNuclearDispatch import IndirectNuclearDispatchOutputs as INDO

class IndirectNuclearTES(NuclearTES): 
    """
    The IndirectNuclearTES class intializes, updates, and runs SSC simulations through PySAM,
    specifically for the SSC tcsmolten_salt module. 
    """
    
    def __init__(self, plant_name="nuclear_mspt_indirect_tes", json_name="model1", **kwargs):
        """ Initializes the SolarTES module
        
        Inputs:
            plant_name (str)         : name of SSC module to run 
            json_name (str)          : name of JSON script with input data for module
            is_dispatch (bool)       : boolean, if True runs Pyomo dispatch optimization
        """
        self.u = SSCHelperMethods.define_unit_registry()
        
        # calculating steam properties from steam table
        steampath = os.path.join( FileMethods.data_dir, "steam_table.csv")
        T, cp, Hp = FileMethods.read_steam_table_file( steampath, self.u )
        self.cp_interp = interp1d( T.m, cp.m, kind='linear' ) # T: in deg K || cp: in kJ/kg*K
        self.hp_interp = interp1d( T.m, Hp.m, kind='linear' ) # T: in dek K || Hp: in kJ/kg
        
        # overriding base class default value unless we get a keyword from higher up
        kwargs['direct'] = kwargs['direct'] if 'direct' in kwargs else False
        super().__init__( plant_name, json_name, **kwargs )
        
        # define specific PySAM module to be called later
        self.PySAM_Module = NuclearMsptIndirectTes
        
        # define specific Dispatch module to be called later
        self.Dispatch_Module = IND
        
        # define the specific Dispatch Outputs module to be called later to create dispatch targets for SSC
        self.Dispatch_Outputs = INDO


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
        
        self.DispatchParameterClass = INDP
        
        dispatch_wrap = self.DispatchParameterClass( unit_registry=self.u, 
                    SSC_dict=self.SSC_dict, PySAM_dict=PySAM_dict,
                    pyomo_horizon=self.pyomo_horizon, 
                    dispatch_time_step=self.dispatch_time_step,
                    interpolants=self.interpolants)
        
        return dispatch_wrap


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
                'q_thermal_log':     'Q_nuc_thermal',    # thermal power from nuclear to HTF 
                'p_cycle_log' :      'P_cycle',          # PC electrical power output (gross)
                'q_dot_rec_inc_log': 'q_dot_nuc_inc',    # Nuclear incident thermal power
                'q_pb_log':          'q_pb',             # PC input energy
                'q_dot_pc_su_log' :  'q_dot_pc_startup', # PC startup thermal power
                'q_dc_tes' :         'q_dc_tes',         # TES discharge thermal power
                'q_ch_tes' :         'q_ch_tes',         # TES charge thermal power
                'm_dot_pc_log' :     'm_dot_pc',         # PC HTF mass flow rate
                'm_dot_rec_log'  :   'm_dot_nuc',        # Nuc mass flow rate
                'T_pc_in_log' :      'T_pc_in',          # PC HTF inlet temperature 
                'T_pc_out_log'   :   'T_pc_out',         # PC HTF outlet temperature
                'T_tes_cold_log':    'T_tes_cold',       # TES cold temperature
                'T_tes_hot_log'  :   'T_tes_hot',        # TES hot temperature
                'T_rec_in_log':      'T_nuc_in',         # Plant inlet temperature
                'T_rec_out_log'  :   'T_nuc_out',        # Plant outlet temperature
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
                   'f_nuc_to_tes_target'  : empty_array.copy(),
                   'is_pc_sb_allowed_in'  : empty_array.copy(),
                   'q_pc_target_su_in'    : empty_array.copy(),
                   'q_pc_target_on_in'    : empty_array.copy(),
                   'q_pc_max_in'          : empty_array.copy()
                   }
            
            
    def create_dispatch_params(self, Plant):
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
        params = NuclearTES.create_dispatch_params(self, Plant)
        
        # specific params for IndirectNuclearDispatch
        params = DW.set_indirect_config_parameters( params )
        
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
        
        Plant.SystemControl.f_nuc_to_tes_target  = dispatch_targets['f_nuc_to_tes_target']
        
        return Plant