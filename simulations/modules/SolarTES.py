#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 13:49:24 2021

@author: gabrielsoto
"""

import sys
sys.path.append('..')
import PySAM.TcsmoltenSalt as TcsmoltenSalt
from modules.NuclearTES import NuclearTES
# from dispatch.NuclearDispatch import NuclearDispatch as ND
# from dispatch.NuclearDispatch import NuclearDispatchParamWrap as NDP

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


    def update_Plant_after_SSC(self):
        """ Update SSC Plant inputs with SSC outputs from previous segment simulation
        
        ** self.run_loop == True
        
        This method uses the SSC end results from the previous simulation segment
        and sets them as the initial conditions for the next SSC segment. As a 
        small note: some outputs are arrays that span the full year, however the
        only relevant parts are the first indeces corresponding to the SSC Horizon.
        All other values are typically 0. 
        """
        self.Plant.SystemControl.rec_op_mode_initial              = self.Plant.Outputs.rec_op_mode_final[self.t_ind-1]
        self.Plant.SystemControl.rec_startup_time_remain_init     = self.Plant.Outputs.rec_startup_time_remain_final[self.t_ind-1]
        self.Plant.SystemControl.rec_startup_energy_remain_init   = self.Plant.Outputs.rec_startup_energy_remain_final[self.t_ind-1]
        self.Plant.SystemControl.T_tank_cold_init                 = self.Plant.Outputs.T_tes_cold[self.t_ind-1]
        self.Plant.SystemControl.T_tank_hot_init                  = self.Plant.Outputs.T_tes_hot[self.t_ind-1]
        self.Plant.ThermalStorage.csp_pt_tes_init_hot_htf_percent = self.Plant.Outputs.hot_tank_htf_percent_final[self.t_ind-1]
        self.Plant.SystemControl.pc_op_mode_initial               = self.Plant.Outputs.pc_op_mode_final[self.t_ind-1]
        self.Plant.SystemControl.pc_startup_energy_remain_initial = self.Plant.Outputs.pc_startup_time_remain_final[self.t_ind-1]
        self.Plant.SystemControl.pc_startup_time_remain_init      = self.Plant.Outputs.pc_startup_energy_remain_final[self.t_ind-1]
        self.Plant.SystemControl.is_field_tracking_init           = self.Plant.Outputs.is_field_tracking_final[self.t_ind-1]
