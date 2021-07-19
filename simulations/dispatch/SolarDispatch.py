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


    def add_objective(self):
        """ Method to add an objective function to the Pyomo Solar Model
        
        This method adds an objective function to the Pyomo Solar Dispatch
        model. Typically, the LORE team defined a nested function and passes it
        into the Pyomo model. 
        """
        def objectiveRule(model):
            """ Maximize profits vs. costs """
            return (
                    sum( model.D[t] * 
                    #obj_profit
                    model.Delta[t]*model.P[t]*0.1*(model.wdot_s[t] - model.wdot_p[t])
                    #obj_cost_cycle_su_hs_sd
                    - (model.Ccsu*model.ycsup[t] + 0.1*model.Cchsp*model.ychsp[t] + model.alpha*model.ycsd[t])
                    #obj_cost_cycle_ramping
                    - (model.C_delta_w*(model.wdot_delta_plus[t]+model.wdot_delta_minus[t])+model.C_v_w*(model.wdot_v_plus[t] + model.wdot_v_minus[t]))
                    #obj_cost_rec_su_hs_sd
                    - (model.Crsu*model.yrsup[t] + model.Crhsp*model.yrhsp[t] + model.alpha*model.yrsd[t])
                    #obj_cost_ops
                    - model.Delta[t]*(model.Cpc*model.wdot[t] + model.Ccsb*model.Qb*model.ycsb[t] + model.Crec*model.xr[t] )
                    for t in model.T) 
                    )
        
        self.model.OBJ = pe.Objective(rule=objectiveRule, sense = pe.maximize)


    def addReceiverStartupConstraints(self):
        """ Method to add solar startup constraints to the Pyomo Solar Model
        
        This method adds constraints pertaining to CSP power tower startup within the 
        Pyomo Solar Dispatch class. Several nested functions are defined. 
        """
        def rec_inventory_rule(model,t):
            """ CSP energy inventory for startup is balanced """
            if t == 1:
                return model.ursu[t] <= model.ursu0 + model.Delta[t]*model.xrsu[t]
            return model.ursu[t] <= model.ursu[t-1] + model.Delta[t]*model.xrsu[t]
        def rec_inv_nonzero_rule(model,t):
            """ CSP energy inventory is positive if starting up """
            return model.ursu[t] <= model.Er * model.yrsu[t]
        def rec_startup_rule(model,t):
            """ CSP resumes normal operation after startup or if previously operating """
            if t == 1:
                return model.yr[t] <= model.ursu[t]/model.Er + model.yr0 + model.yrsb0
            return model.yr[t] <= model.ursu[t]/model.Er + model.yr[t-1] + model.yrsb[t-1]
        def rec_su_persist_rule(model,t):
            """ CSP cannot startup if operating in previous timestep """
            if t == 1: 
                return model.yrsu[t] + model.yr0 <= 1
            return model.yrsu[t] +  model.yr[t-1] <= 1
        def ramp_limit_rule(model,t):
            """ CSP must abide by ramp rate limits """
            return model.xrsu[t] <= model.Qru*model.yrsu[t]
        def nontrivial_solar_rule(model,t):
            """ CSP trivial power prevents startup """
            return model.yrsu[t] <= model.Qin[t]
        
        self.model.rec_inventory_con = pe.Constraint(self.model.T,rule=rec_inventory_rule)
        self.model.rec_inv_nonzero_con = pe.Constraint(self.model.T,rule=rec_inv_nonzero_rule)
        self.model.rec_startup_con = pe.Constraint(self.model.T,rule=rec_startup_rule)
        self.model.rec_su_persist_con = pe.Constraint(self.model.T,rule=rec_su_persist_rule)
        self.model.ramp_limit_con = pe.Constraint(self.model.T,rule=ramp_limit_rule)
        self.model.nontrivial_solar_con = pe.Constraint(self.model.T,rule=nontrivial_solar_rule)


    def addReceiverSupplyAndDemandConstraints(self):
        """ Method to add CSP supply and demand constraints to the Pyomo Solar Model
        
        This method adds constraints pertaining to solar supply and demand energy
        constraints. 
        """
        def rec_production_rule(model,t):
            """ Upper bound on thermal energy produced by NP """
            return model.xr[t] + model.xrsu[t] + model.Qrsd*model.yrsd[t] <= model.Qin[t]
        def rec_generation_rule(model,t):
            """ Thermal energy production by NP only when operating """
            return model.xr[t] <= model.Qin[t] * model.yr[t]
        def min_generation_rule(model,t):
            """ Lower bound on thermal energy produced by NP """
            return model.xr[t] >= model.Qrl * model.yr[t]
        def rec_gen_persist_rule(model,t):
            """ NP not able to operate if no thermal power (adapted from CSP) """
            return model.yr[t] <= model.Qin[t]/model.Qrl
        
        self.model.rec_production_con = pe.Constraint(self.model.T,rule=rec_production_rule)
        self.model.rec_generation_con = pe.Constraint(self.model.T,rule=rec_generation_rule)
        self.model.min_generation_con = pe.Constraint(self.model.T,rule=min_generation_rule)
        self.model.rec_gen_persist_con = pe.Constraint(self.model.T,rule=rec_gen_persist_rule)
        

    def addReceiverNodeLogicConstraints(self):
        """ Method to add solar logic constraints to the Pyomo Solar Dispatch Model
        
        This method adds constraints pertaining to tracking binary variable
        logic for the CSP. Essentially, to make sure the correct modes
        are activated when they are allowed. 
        """
        def rec_su_sb_persist_rule(model,t):
            """ CSP startup and standby cannot coincide """
            return model.yrsu[t] + model.yrsb[t] <= 1
        def rec_sb_persist_rule(model,t):
            """ CSP standby and operation cannot coincide """
            return model.yr[t] + model.yrsb[t] <= 1
        def rsb_persist_rule(model,t):
            """ CSP standby must follow operating or another standby timestep """
            if t == 1:
                return model.yrsb[t] <= (model.yr0 + model.yrsb0) 
            return model.yrsb[t] <= model.yr[t-1] + model.yrsb[t-1]
        def rec_su_pen_rule(model,t):
            """ CSP cold startup penalty """
            if t == 1:
                return model.yrsup[t] >= model.yrsu[t] - model.yrsu0 
            return model.yrsup[t] >= model.yrsu[t] - model.yrsu[t-1]
        def rec_hs_pen_rule(model,t):
            """ CSP hot startup penalty: if we went from standby to operating """
            if t == 1:
                return model.yrhsp[t] >= model.yr[t] - (1 - model.yrsb0)
            return model.yrhsp[t] >= model.yr[t] - (1 - model.yrsb[t-1])
        def rec_shutdown_rule(model,t):
            """ CSP shutdown protocol """
            current_Delta = model.Delta[t]
            # structure of inequality is lb <= model.param <= ub with strict=False by default
            if self.eval_ineq(1,current_Delta) and t == 1: #not strict
                return 0 >= model.yr0 - model.yr[t] +  model.yrsb0 - model.yrsb[t]
            elif self.eval_ineq(1,current_Delta) and t > 1: # not strict
                return model.yrsd[t-1] >= model.yr[t-1] - model.yr[t] + model.yrsb[t-1] - model.yrsb[t]
            elif self.eval_ineq(current_Delta,1,strict=True) and t == 1:
                return model.yrsd[t] >= model.yr0  - model.yr[t] + model.yrsb0 - model.yrsb[t]
            # only case remaining: Delta[t]<1, t>1
            return model.yrsd[t] >= model.yr[t-1] - model.yr[t] + model.yrsb[t-1] - model.yrsb[t]
        
        self.model.rec_su_sb_persist_con = pe.Constraint(self.model.T,rule=rec_su_sb_persist_rule)
        self.model.rec_sb_persist_con = pe.Constraint(self.model.T,rule=rec_sb_persist_rule)
        self.model.rsb_persist_con = pe.Constraint(self.model.T,rule=rsb_persist_rule)
        self.model.rec_su_pen_con = pe.Constraint(self.model.T,rule=rec_su_pen_rule)
        self.model.rec_hs_pen_con = pe.Constraint(self.model.T,rule=rec_hs_pen_rule)
        self.model.rec_shutdown_con = pe.Constraint(self.model.T,rule=rec_shutdown_rule)