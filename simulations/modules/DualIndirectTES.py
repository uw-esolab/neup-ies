#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 16:28:05 2022

@author: gabrielsoto
"""

import sys, copy
sys.path.append('..')
import numpy as np
import PySAM.NuclearMsptIndirectTes as NuclearMsptIndirectTes
from modules.GenericSSCModule import GenericSSCModule
from modules.IndirectNuclearTES import IndirectNuclearTES
from modules.SolarTES import SolarTES
from dispatch.DualIndirectDispatch import DualIndirectDispatch as DI
from dispatch.DualIndirectDispatch import DualIndirectDispatchParamWrap as DIP
from dispatch.DualIndirectDispatch import DualIndirectDispatchOutputs as DIO
from dispatch.SolarDispatch import SolarDispatchParamWrap as SDP

class DualIndirectTES(SolarTES): 
    """
    The DualPlantTES class intializes, updates, and runs SSC simulations through PySAM,
    specifically for the SSC nuclear_mspt_tes module. 
    
    This is meant to simulate the intermediate cycle where the Nuclear and Solar 
    plants are both directly connected to the storage tank - power cycle loop. 
    That is, the Nuclear and Solar Power Tower heat parallel mass flows
    of molten salt via respective heat exchangers. Each molten salt mass flow
    then routes directly to the hot storage tank where it can be dispatched
    out to the power cycle. 
    """

    def __init__(self, plant_name="nuclear_mspt_tes", json_name="model2", **specs):
        """ Initializes the DualPlantTES module
        
        Inputs:
            plant_name (str)         : name of SSC module to run 
            json_name (str)          : name of JSON script with input data for module
            is_dispatch (bool)       : boolean, if True runs Pyomo dispatch optimization
        """
        
        # initialize Solar+Nuclear+Generic module, csv data arrays should be saved here
        SolarTES.__init__( self, plant_name, json_name, **specs )
        
        # define specific PySAM module to be called later
        self.PySAM_Module = NuclearMsptIndirectTes

        # define specific Dispatch module to be called later
        self.Dispatch_Module = DI
        
        self.Dispatch_Outputs = DIO


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
        
        self.DispatchParameterClass = DIP
        
        dispatch_wrap = self.DispatchParameterClass( unit_registry=self.u, 
                    SSC_dict=self.SSC_dict, PySAM_dict=PySAM_dict,
                    pyomo_horizon=self.pyomo_horizon, 
                    dispatch_time_step=self.dispatch_time_step)
        
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
        
        # extract array from full run of Plant for solar
        assert hasattr(Plant.Outputs, "Q_thermal"), "Q_thermal was not found in the outputs of Plant."
        self.Q_rec_guess = Plant.Outputs.Q_thermal

        # extract array from full run of Plant
        assert hasattr(Plant.Outputs, "Q_nuc_thermal"), "Q_nuc_thermal was not found in the outputs of Plant."
        self.Q_nuc_guess = Plant.Outputs.Q_nuc_thermal
        
        # set up copy of SSC dict
        updated_SSC_dict = copy.deepcopy(self.SSC_dict)
        updated_SSC_dict['Q_thermal']     = self.Q_rec_guess[self.slice_pyo_firstH]
        updated_SSC_dict['Q_nuc_thermal'] = self.Q_nuc_guess[self.slice_pyo_firstH]
        
        # these are NuclearTES-specific setters
        params = DW.set_nuclear_parameters( params )
        params = DW.set_time_series_nuclear_parameters( params, updated_SSC_dict )
        
        # specific params for IndirectNuclearDispatch
        params = DW.set_indirect_config_parameters( params )
        
        # these are SolarTES-specific setters
        params = SDP.set_solar_parameters( DW, params )
        params = SDP.set_time_series_solar_parameters( DW, params, updated_SSC_dict )
        
        # this sets the initial set for the NuclearTES
        params = DW.set_initial_state( params )
        
        return params


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
        
        # field and receiver initial conditions
        self.Plant.SystemControl.rec_op_mode_initial              = self.Plant.Outputs.rec_op_mode_final
        self.Plant.SystemControl.rec_startup_time_remain_init     = self.Plant.Outputs.rec_startup_time_remain_final
        self.Plant.SystemControl.rec_startup_energy_remain_init   = self.Plant.Outputs.rec_startup_energy_remain_final
        self.Plant.SystemControl.is_field_tracking_init           = self.Plant.Outputs.is_field_tracking_final
        
        # nuclear initial conditions
        self.Plant.SystemControl.nuc_op_mode_initial              = self.Plant.Outputs.nuc_op_mode_final
        self.Plant.SystemControl.nuc_startup_time_remain_init     = self.Plant.Outputs.nuc_startup_time_remain_final
        self.Plant.SystemControl.nuc_startup_energy_remain_init   = self.Plant.Outputs.nuc_startup_energy_remain_final
        
        # TES initial conditions
        self.Plant.SystemControl.T_tank_cold_init                 = self.Plant.Outputs.T_tes_cold[ssc_slice][-1]
        self.Plant.SystemControl.T_tank_hot_init                  = self.Plant.Outputs.T_tes_hot[ssc_slice][-1]
        self.Plant.ThermalStorage.csp_pt_tes_init_hot_htf_percent = self.Plant.Outputs.hot_tank_htf_percent_final
        
        # PC initial conditions
        self.Plant.SystemControl.pc_op_mode_initial               = self.Plant.Outputs.pc_op_mode_final
        self.Plant.SystemControl.pc_startup_energy_remain_initial = self.Plant.Outputs.pc_startup_energy_remain_final
        self.Plant.SystemControl.pc_startup_time_remain_init      = self.Plant.Outputs.pc_startup_time_remain_final


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
        updated_SSC_dict['rec_op_mode_initial']              = Plant.Outputs.rec_op_mode_final
        updated_SSC_dict['rec_startup_time_remain_init']     = Plant.Outputs.rec_startup_time_remain_final
        updated_SSC_dict['rec_startup_energy_remain_init']   = Plant.Outputs.rec_startup_energy_remain_final
        updated_SSC_dict['nuc_op_mode_initial']              = Plant.Outputs.nuc_op_mode_final
        updated_SSC_dict['nuc_startup_time_remain_init']     = Plant.Outputs.nuc_startup_time_remain_final
        updated_SSC_dict['nuc_startup_energy_remain_init']   = Plant.Outputs.nuc_startup_energy_remain_final
        updated_SSC_dict['T_tank_cold_init']                 = Plant.Outputs.T_tes_cold[ssc_slice][-1]
        updated_SSC_dict['T_tank_hot_init']                  = Plant.Outputs.T_tes_hot[ssc_slice][-1]
        updated_SSC_dict['csp.pt.tes.init_hot_htf_percent']  = Plant.Outputs.hot_tank_htf_percent_final
        updated_SSC_dict['pc_op_mode_initial']               = Plant.Outputs.pc_op_mode_final
        updated_SSC_dict['pc_startup_time_remain_init']      = Plant.Outputs.pc_startup_time_remain_final
        updated_SSC_dict['pc_startup_energy_remain_initial'] = Plant.Outputs.pc_startup_energy_remain_final
        updated_SSC_dict['is_field_tracking_init']           = Plant.Outputs.is_field_tracking_final
        
        # these are specific to the initial states
        updated_SSC_dict['wdot0']     = Plant.Outputs.P_cycle[ssc_slice][-1]
        
        # extract time series from a previous SSC run
        updated_SSC_dict['Q_thermal']     = self.Q_rec_guess[self.slice_pyo_currentH]
        updated_SSC_dict['Q_nuc_thermal'] = self.Q_nuc_guess[self.slice_pyo_currentH]
        
        # TODO: removing w_dot_s_prev references in all of Dispatch for now, might need to revisit later
        # updated_SSC_dict['wdot_s_prev'] = 0 #np.array([pe.value(dm.model.wdot_s_prev[t]) for t in dm.model.T])[-1]
        
        DW = self.dispatch_wrap
        
        # set up Plant outputs dictionary
        plant_dict = Plant.Outputs.export()
        
        # updating the initial state and time series Nuclear params
        params = DW.set_time_indexed_parameters( params, self.df_array, self.ud_array, self.slice_pyo_currentH )
        params = DW.set_initial_state( params, updated_SSC_dict, Plant, self.t_ind )
        params = DW.set_time_series_nuclear_parameters( params, updated_SSC_dict )
        params = SDP.set_time_series_solar_parameters( DW, params, updated_SSC_dict )
        
        return params