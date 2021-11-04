#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:03:51 2021

@author: gabrielsoto
"""

import sys, copy
sys.path.append('..')
import numpy as np
import PySAM.NuclearMsptTes as NuclearMsptTes
from modules.GenericSSCModule import GenericSSCModule
from modules.NuclearTES import NuclearTES
from modules.SolarTES import SolarTES

class DualPlantTES(SolarTES): 
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

    def __init__(self, plant_name="nuclear_mspt_tes", json_name="model2", is_dispatch=False):
        """ Initializes the DualPlantTES module
        
        Inputs:
            plant_name (str)         : name of SSC module to run 
            json_name (str)          : name of JSON script with input data for module
            is_dispatch (bool)       : boolean, if True runs Pyomo dispatch optimization
        """
        
        # initialize Solar+Nuclear+Generic module, csv data arrays should be saved here
        SolarTES.__init__( self, plant_name, json_name, is_dispatch )
        
        # define specific PySAM module to be called later
        self.PySAM_Module = NuclearMsptTes