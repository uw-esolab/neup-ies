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
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
import copy

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