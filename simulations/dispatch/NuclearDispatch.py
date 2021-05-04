# -*- coding: utf-8 -*-
"""
Pyomo real-time dispatch model


Model authors:
* Mike Wagner - UW-Madison
* Bill Hamilton - NREL 
* John Cox - Colorado School of Mines
* Alex Zolan - NREL

Pyomo code by Alex Zolan
Modified by Gabriel Soto
"""
from dispatch.GeneralDispatch import GeneralDispatch
from dispatch.GeneralDispatch import GeneralDispatchParamWrap
import pyomo.environ as pe
import numpy as np
import pint
u = pint.UnitRegistry()

class NuclearDispatch(GeneralDispatch):
    
    def __init__(self, params):
        # initialize Generic module, csv data arrays should be saved here
        GeneralDispatch.__init__( self, params )

# =============================================================================
# Dispatch Wrapper
# =============================================================================

class NuclearDispatchParamWrap(GeneralDispatchParamWrap):
    
    def __init__(self, SSC_dict=None, PySAM_dict=None, pyomo_horizon=48, 
                   dispatch_time_step=1):
        
        GeneralDispatchParamWrap.__init__( self, SSC_dict, PySAM_dict, 
                            pyomo_horizon, dispatch_time_step )


    def set_fixed_cost_parameters(self, param_dict):
    
            # set up costs from parent class
            param_dict = GeneralDispatchParamWrap.set_fixed_cost_parameters( param_dict )
            
            return param_dict