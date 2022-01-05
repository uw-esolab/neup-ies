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
import pyomo.environ as pe
from dispatch.NuclearDispatch import NuclearDispatch
from dispatch.SolarDispatch   import SolarDispatch
from dispatch.GeneralDispatch import GeneralDispatchParamWrap
import numpy as np
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
import os, copy

class DualPlantDispatch(NuclearDispatch):
    """
    The DualPlantDispatch class is meant to set up and run Dispatch
    optimization as a mixed integer linear program problem using Pyomo,
    specifically for the NuclearMsptTES NE2+SSC module.
    """

    def __init__(self, params, unitRegistry):
        """ Initializes the DualPlantDispatch module
        
        The instantiation of this class receives a parameter dictionary from
        the NE2 module (created using the DualPlantDispatchWrapper class). It calls
        on the GeneralDispatch __init__ to create the model. The NuclearDispatch first
        creates an empty Concrete Model from Pyomo, then generates Parameters
        from the parameter dictionary, Variables, Objectives and Constraints.
        
        Inputs:
            params (dict)                : dictionary of Pyomo dispatch parameters
            unitRegistry (pint.registry) : unique unit Pint unit registry
        """
        
        # initialize Nuclear module, csv data arrays should be saved here
        NuclearDispatch.__init__( self, params, unitRegistry )


    def generate_params(self, params):
        """ Method to generate parameters within Pyomo DualPlant Model
        
        This method reads in a dictionary of pyomo parameters and uses these
        inputs to initialize parameters for the Pyomo Concrete Model. This method
        sets up parameters particularly for the Power Cycle. It also defines
        some lambda functions that helps convert Pint units to Pyomo units. It
        first instantiates PowerCycle parameters through GeneralDispatch, then
        instantiates Nuclear parameters.
        
        Note: initial conditions are defined for the time period immediately 
        preceding the start of this new Pyomo time segment. 
        
        Inputs:
            params (dict)  : dictionary of Pyomo dispatch parameters
        """
        
        # generating GeneralDispatch parameters first (PowerCycle, etc.)
        NuclearDispatch.generate_params(self, params)
        SolarDispatch.generate_params(self, params, skip_parent=True)


    def generate_variables(self):
        """ Method to generate parameters within Pyomo DualPlant Model
        
        This method instantiates variables for the Pyomo Concrete Model, with
        domains. Does not need initial guesses here, they are defined in the 
        parameters. We first define continuous and binary variables for the 
        Power Cycle through GeneralDispatch, then declare nuclear variables.
        """
        
        # generating GeneralDispatch variables first (PowerCycle, etc.)
        NuclearDispatch.generate_variables(self)
        SolarDispatch.generate_params(self, skip_parent=True)