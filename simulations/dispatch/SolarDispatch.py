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
from dispatch.GeneralDispatch import GeneralDispatch
from dispatch.GeneralDispatch import GeneralDispatchParamWrap
import numpy as np
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
import os, copy

class SolarDispatch(GeneralDispatch):
    """
    The SolarDispatch class is meant to set up and run Dispatch
    optimization as a mixed integer linear program problem using Pyomo,
    specifically for the NuclearTES NE2+SSC module.
    """
    
    def __init__(self, params, unitRegistry):
        """ Initializes the SolarDispatch module
        
        The instantiation of this class receives a parameter dictionary from
        the NE2 module (created using the SolarDispatchWrapper class). It calls
        on the GeneralDispatch __init__ to create the model. The GeneralDispatcher first
        creates an empty Concrete Model from Pyomo, then generates Parameters
        from the parameter dictionary, Variables, Objectives and Constraints.
        
        Inputs:
            params (dict)                : dictionary of Pyomo dispatch parameters
            unitRegistry (pint.registry) : unique unit Pint unit registry
        """
        
        # initialize Generic module, csv data arrays should be saved here
        GeneralDispatch.__init__( self, params, unitRegistry )


    def generate_params(self, params):
        """ Method to generate parameters within Pyomo Solar Model
        
        This method reads in a dictionary of pyomo parameters and uses these
        inputs to initialize parameters for the Pyomo Concrete Model. This method
        sets up parameters particularly for the Power Cycle. It also defines
        some lambda functions that helps convert Pint units to Pyomo units. It
        first instantiates PowerCycle parameters through GeneralDispatch, then
        instantiates Solar parameters.
        
        Note: initial conditions are defined for the time period immediately 
        preceding the start of this new Pyomo time segment. 
        
        Inputs:
            params (dict)  : dictionary of Pyomo dispatch parameters
        """
        
        # generating GeneralDispatch parameters first (PowerCycle, etc.)
        GeneralDispatch.generate_params(self, params)
        
        # lambdas to convert units and data to proper syntax
        gd = self.gd
        gu = self.gu 
        
        ### Cost Parameters ### 
        self.model.Crec = pe.Param(mutable=True, initialize=gd("Crec"), units=gu("Crec"))              #C^{rec}: Operating cost of heliostat field and receiver [\$/kWt$\cdot$h]
        self.model.Crsu = pe.Param(mutable=True, initialize=gd("Crsu"), units=gu("Crsu"))              #C^{rsu}: Penalty for receiver cold start-up [\$/start]
        self.model.Crhsp = pe.Param(mutable=True, initialize=gd("Crhsp"), units=gu("Crhsp"))           #C^{rhsp}: Penalty for receiver hot start-up [\$/start]
        
        ### CSP Parameters ###
        #TODO: Eu is reused here
        self.model.deltal = pe.Param(mutable=True, initialize=gd("deltal"), units=gu("deltal"))    #\delta^l: Minimum time to start the receiver [hr]
        self.model.Ehs = pe.Param(mutable=True, initialize=gd("Ehs"), units=gu("Ehs"))             #E^{hs}: Heliostat field startup or shut down parasitic loss [kWe$\cdot$h]
        self.model.Er = pe.Param(mutable=True, initialize=gd("Er"), units=gu("Er"))                #E^r: Required energy expended to start receiver [kWt$\cdot$h]
        self.model.Eu = pe.Param(mutable=True, initialize=gd("Eu"), units=gu("Eu"))                #E^u: Thermal energy storage capacity [kWt$\cdot$h]
        self.model.Lr = pe.Param(mutable=True, initialize=gd("Lr"), units=gu("Lr"))                #L^r: Receiver pumping power per unit power produced [kWe/kWt]
        self.model.Qrl = pe.Param(mutable=True, initialize=gd("Qrl"), units=gu("Qrl"))             #Q^{rl}: Minimum operational thermal power delivered by receiver [kWt$\cdot$h]
        self.model.Qrsb = pe.Param(mutable=True, initialize=gd("Qrsb"), units=gu("Qrsb"))          #Q^{rsb}: Required thermal power for receiver standby [kWt$\cdot$h]
        self.model.Qrsd = pe.Param(mutable=True, initialize=gd("Qrsd"), units=gu("Qrsd"))          #Q^{rsd}: Required thermal power for receiver shut down [kWt$\cdot$h] 
        self.model.Qru = pe.Param(mutable=True, initialize=gd("Qru"), units=gu("Qru"))             #Q^{ru}: Allowable power per period for receiver start-up [kWt$\cdot$h]
        self.model.Wh = pe.Param(mutable=True, initialize=gd("Wh"), units=gu("Wh"))                #W^h: Heliostat field tracking parasitic loss [kWe]
        self.model.Wht = pe.Param(mutable=True, initialize=gd("Wht"), units=gu("Wht"))             #W^{ht}: Tower piping heat trace parasitic loss [kWe]

        ### Time series CSP Parameters ###
        self.model.delta_rs = pe.Param(self.model.T, mutable=True, initialize=gd("delta_rs"), units=gu("delta_rs"))     #\delta^{rs}_{t}: Estimated fraction of period $t$ required for receiver start-up [-]
        self.model.Qin = pe.Param(self.model.T, mutable=True, initialize=gd("Qin"), units=gu("Qin"))                    #Q^{in}_{t}: Available thermal power generated by the CSP heliostat field in period $t$ [kWt]

        ### Initial Condition Parameters ###
        self.model.s0 = pe.Param(mutable=True, initialize=gd("s0"), units=gu("s0"))           #s_0: Initial TES reserve quantity  [kWt$\cdot$h]
        self.model.ursu0 = pe.Param(mutable=True, initialize=gd("ursu0"), units=gu("ursu0"))  #u^{rsu}_0: Initial receiver start-up energy inventory [kWt$\cdot$h]
        self.model.yr0 = pe.Param(mutable=True, initialize=gd("yr0"), units=gu("yr0"))        #y^r_0: 1 if receiver is generating ``usable'' thermal power initially, 0 otherwise  [az] this is new.
        self.model.yrsb0 = pe.Param(mutable=True, initialize=gd("yrsb0"), units=gu("yrsb0"))  #y^{rsb}_0: 1 if receiver is in standby mode initially, 0 otherwise [az] this is new.
        self.model.yrsu0 = pe.Param(mutable=True, initialize=gd("yrsu0"), units=gu("yrsu0"))  #y^{rsu}_0: 1 if receiver is in starting up initially, 0 otherwise    [az] this is new.
        

    def generate_variables(self):
        """ Method to generate parameters within Pyomo Solar Model
        
        This method instantiates variables for the Pyomo Concrete Model, with
        domains. Does not need initial guesses here, they are defined in the 
        parameters. We first define continuous and binary variables for the 
        Power Cycle through GeneralDispatch, then declare nuclear variables.
        """
        
        # generating GeneralDispatch variables first (PowerCycle, etc.)
        GeneralDispatch.generate_variables(self)
        
        ### Decision Variables ###
        #------- Variables ---------
        self.model.s = pe.Var(self.model.T, domain=pe.NonNegativeReals, bounds = (0,self.model.Eu))    #s: TES reserve quantity at period $t$  [kWt$\cdot$h]
        self.model.ursu = pe.Var(self.model.T, domain=pe.NonNegativeReals)     #u^{rsu}: Receiver start-up energy inventory at period $t$ [kWt$\cdot$h]
        self.model.xr = pe.Var(self.model.T, domain=pe.NonNegativeReals)       #x^r: Thermal power delivered by the receiver at period $t$ [kWt]
        self.model.xrsu = pe.Var(self.model.T, domain=pe.NonNegativeReals)     #x^{rsu}: Receiver start-up power consumption at period $t$ [kWt]
        
        #------- Binary Variables ---------
        self.model.yr = pe.Var(self.model.T, domain=pe.Binary)        #y^r: 1 if receiver is generating ``usable'' thermal power at period $t$; 0 otherwise
        self.model.yrhsp = pe.Var(self.model.T, domain=pe.Binary)	  #y^{rhsp}: 1 if receiver hot start-up penalty is incurred at period $t$ (from standby); 0 otherwise
        self.model.yrsb = pe.Var(self.model.T, domain=pe.Binary)	  #y^{rsb}: 1 if receiver is in standby mode at period $t$; 0 otherwise
        self.model.yrsd = pe.Var(self.model.T, domain=pe.Binary)	  #y^{rsd}: 1 if receiver is shut down at period $t$; 0 otherwise
        self.model.yrsu = pe.Var(self.model.T, domain=pe.Binary)      #y^{rsu}: 1 if receiver is starting up at period $t$; 0 otherwise
        self.model.yrsup = pe.Var(self.model.T, domain=pe.Binary)     #y^{rsup}: 1 if receiver cold start-up penalty is incurred at period $t$ (from off); 0 otherwise

