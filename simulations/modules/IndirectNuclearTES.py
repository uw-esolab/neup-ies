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
        self.cp_interp = interp1d( T, cp, kind='linear' )
        self.hp_interp = interp1d( T, Hp, kind='linear' )
        
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
                    dispatch_time_step=self.dispatch_time_step)
        
        return dispatch_wrap


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