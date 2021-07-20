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


    def addTESEnergyBalanceConstraints(self):
        """ Method to add TES constraints to the Pyomo Solar Model
        
        This method adds constraints pertaining to TES energy balance from charging
        with thermal power and discharging to the power cycle. 
        """
        def tes_balance_rule(model, t):
            """ Balance of energy to and from TES """
            if t == 1:
                return model.s[t] - model.s0 == model.Delta[t] * (model.xr[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.x[t] + model.Qrsb*model.yrsb[t]))
            return model.s[t] - model.s[t-1] == model.Delta[t] * (model.xr[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.x[t] + model.Qrsb*model.yrsb[t]))
        def tes_upper_rule(model, t):
            """ Upper bound to TES charge state """
            return model.s[t] <= model.Eu
        def tes_start_up_rule(model, t):
            """ Ensuring sufficient TES charge level to startup NP """
            if t == 1:
                return model.s0 >= model.Delta[t]*model.delta_rs[t]*( (model.Qu + model.Qb)*( -3 + model.yrsu[t] + model.y0 + model.y[t] + model.ycsb0 + model.ycsb[t] ) + model.x[t] + model.Qb*model.ycsb[t] )
            return model.s[t-1] >= model.Delta[t]*model.delta_rs[t]*( (model.Qu + model.Qb)*( -3 + model.yrsu[t] + model.y[t-1] + model.y[t] + model.ycsb[t-1] + model.ycsb[t] ) + model.x[t] + model.Qb*model.ycsb[t] )
        def maintain_tes_rule(model):
            """ Final state of TES has to be less than or equal to start """
            return model.s[model.num_periods] <= model.s0
        
        self.model.tes_balance_con = pe.Constraint(self.model.T,rule=tes_balance_rule)
        self.model.tes_upper_con = pe.Constraint(self.model.T,rule=tes_upper_rule)
        self.model.tes_start_up_con = pe.Constraint(self.model.T,rule=tes_start_up_rule)
        self.model.maintain_tes_con = pe.Constraint(rule=maintain_tes_rule)


    def addPiecewiseLinearEfficiencyConstraints(self):
        """ Method to add efficiency constraints to the Pyomo Solar Model
        
        This method adds constraints pertaining to efficiency constraints defined
        as a piecewise linear approximation. Also referred to as Cycle supply and 
        demand constraints. In the SolarDispatch, we add an extra balance of power
        with respect to energy storage and power produced from the CSP plant. 
        
        TODO: This should be revisited when adding MED!!
        """
        def grid_sun_rule(model, t):
            """ Balance of power flow, i.e. sold vs purchased """
            return (
                    model.wdot_s[t] - model.wdot_p[t] == (1-model.etac[t])*model.wdot[t]
                		- model.Lr*(model.xr[t] + model.xrsu[t] + model.Qrl*model.yrsb[t])
                		- model.Lc*model.x[t] 
                        - model.Wh*model.yr[t] - model.Wb*model.ycsb[t] - model.Wht*(model.yrsb[t]+model.yrsu[t])		#Is Wrsb energy [kWh] or power [kW]?  [az] Wrsb = Wht in the math?
                		- (model.Ehs/model.Delta[t])*(model.yrsu[t] + model.yrsb[t] + model.yrsd[t])
            )
        
        # call the parent version of this method
        GeneralDispatch.addPiecewiseLinearEfficiencyConstraints(self)
        
        # additional constraints
        self.model.grid_sun_con = pe.Constraint(self.model.T,rule=grid_sun_rule)
        

    def generate_constraints(self):
        """ Method to add ALL constraints to the Pyomo Solar Model
        
        This method calls the previously defined constraint methods to instantiate
        them and add to the existing model. This method first calls the GeneralDispatch
        version to set PowerCycle constraints, then calls nuclear constraint methods
        to add them to the model. 
        """
        
        # generating GeneralDispatch constraints first (PowerCycle, etc.)
        GeneralDispatch.generate_constraints(self)
        
        self.addReceiverStartupConstraints()
        self.addReceiverSupplyAndDemandConstraints()
        self.addReceiverNodeLogicConstraints()
        self.addTESEnergyBalanceConstraints()


# =============================================================================
# Dispatch Wrapper
# =============================================================================

class SolarDispatchParamWrap(GeneralDispatchParamWrap):
    """
    The SolarDispatchParamWrap class is meant to be the staging area for the 
    creation of Parameters ONLY for the SolarDispatch class. It communicates 
    with the NE2 modules, receiving SSC and PySAM input dictionaries to calculate 
    both static parameters used for every simulation segment AND initial conditions 
    that can be updated.
    """
    
    def __init__(self, unit_registry, SSC_dict=None, PySAM_dict=None, pyomo_horizon=48, 
                   dispatch_time_step=1):
        """ Initializes the SolarDispatchParamWrap module
        
        Inputs:
            unitRegistry (pint.registry)   : unique unit Pint unit registry
            SSC_dict (dict)                : dictionary of SSC inputs needed to run modules
            PySAM_dict (dict)              : dictionary of PySAM inputs + file names
            pyomo_horizon (int Quant)      : length of Pyomo simulation segment (hours)
            dispatch_time_step (int Quant) : length of each Pyomo time step (hours)
        """
        
        GeneralDispatchParamWrap.__init__( self, unit_registry, SSC_dict, PySAM_dict, 
                            pyomo_horizon, dispatch_time_step )


    def set_design(self):
        """ Method to calculate and save design point values of Plant operation
        
        This method extracts values and calculates for design point parameters 
        of our Plant (e.g., nuclear thermal power output, power cycle efficiency,
        inlet and outlet temperatures, etc.). 
        """
        
        u = self.u
        
        GeneralDispatchParamWrap.set_design(self)

        # CSP parameters
        self.q_rec_design = self.SSC_dict['q_dot_rec_inc'] * u.MW      # CSP design thermal power
        
        # specific heat values at design point
        T_htf  = 0.5*(self.T_htf_hot + self.T_htf_cold)
        cp_des = SSCHelperMethods.get_cp_htf(self.u, T_htf, self.SSC_dict['rec_htf'] )
        cp_des = cp_des.to('J/g/kelvin')       
        
        # mass flow rate
        dm_des = self.q_pb_design / (cp_des * (self.T_htf_hot - self.T_htf_cold) )  
        self.dm_pb_design = dm_des.to('kg/s')                               # power block design mass flow rate
        
        # TES design point
        e_tes_design = self.q_pb_design * self.SSC_dict['tshours']*u.hr  
        m_tes_des = e_tes_design / cp_des / (self.T_htf_hot - self.T_htf_cold)     
        self.e_tes_design = e_tes_design.to('kWh') # TES storage capacity (kWht)
        self.m_tes_design = m_tes_des.to('kg')     # TES active storage mass (kg)


    def set_fixed_cost_parameters(self, param_dict):
        """ Method to set fixed costs of the Plant
        
        This method calculates some fixed costs for the Plant operations, startup,
        standby, etc. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        # grabbing unit registry set up in GeneralDispatch
        u = self.u
    
        # set up costs from parent class
        param_dict = GeneralDispatchParamWrap.set_fixed_cost_parameters( self, param_dict )
        
        # TODO: old values from LORE files
        C_rec  = self.PySAM_dict['rec_op_cost'] * u.USD / u.MWh #Q_ratio * 0.002  * u.USD/u.kWh        
        C_rsu  = self.PySAM_dict['rec_cold_su'] * u.USD
        C_rhsp = self.PySAM_dict['rec_hot_su'] * u.USD

        ### Cost Parameters ###
        param_dict['Crec']   = C_rec.to('USD/kWh')  #C^{rec}: Operating cost of heliostat field and receiver [\$/kWt$\cdot$h]
        param_dict['Crsu']   = C_rsu.to('USD')      #C^{rsu}: Penalty for receiver cold start-up [\$/start]
        param_dict['Crhsp']  = C_rhsp.to('USD')     #C^{rhsp}: Penalty for receiver hot start-up [\$/start]

        return param_dict


    def set_solar_parameters(self, param_dict):
        """ Method to set parameters specific to the Solar Plant for Dispatch optimization
        
        This method calculates some parameters specific to the SolarTES plant
        which are meant to be fixed throughout the simulation. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        # grabbing unit registry set up in GeneralDispatch
        u = self.u 
        
        time_fix = 1*u.hr                  # TODO: we're missing a time term to fix units
        dw_rec_pump             = self.PySAM_dict['dw_rec_pump']*u.MW   # TODO: Pumping parasitic at design point reciever mass flow rate (MWe)
        tower_piping_ht_loss    = self.PySAM_dict['tower_piping_ht_loss']*u.kW   # TODO: Tower piping heat trace full-load parasitic load (kWe) 
        q_rec_standby_fraction  = self.PySAM_dict['q_rec_standby_frac']        # TODO: Nuclear standby energy consumption (fraction of design point thermal power)
        q_rec_shutdown_fraction = self.PySAM_dict['q_rec_shutdown_frac']       # TODO: Nuclear shutdown energy consumption (fraction of design point thermal power)
        
        self.deltal = self.SSC_dict['rec_su_delay']*u.hr
        self.Ehs    = self.SSC_dict['p_start']*u.kWh
        self.Er     = self.SSC_dict['rec_qf_delay'] * self.q_rec_design * time_fix
        self.Eu     = self.SSC_dict['tshours']*u.hr * self.q_pb_design
        self.Lr     = dw_rec_pump / self.q_rec_design
        self.Qrl    = self.SSC_dict['f_rec_min'] * self.q_rec_design * time_fix
        self.Qrsb   = q_rec_standby_fraction  * self.q_rec_design * time_fix
        self.Qrsd   = q_rec_shutdown_fraction * self.q_rec_design * time_fix
        self.Qru    = self.Er / self.deltarl  
        self.Wh     = self.SSC_dict['p_track']*u.kW
        self.Wnht   = tower_piping_ht_loss
        
        ### CSP Field and Receiver Parameters ###
        param_dict['deltal'] = self.deltal.to('hr')    #\delta^l: Minimum time to start the receiver [hr]
        param_dict['Ehs']    = self.Ehs.to('kWh')      #E^{hs}: Heliostat field startup or shut down parasitic loss [kWe$\cdot$h]
        param_dict['Er']     = self.Er.to('kWh')       #E^r: Required energy expended to start receiver [kWt$\cdot$h]
        param_dict['Eu']     = self.Eu.to('kWh')       #E^u: Thermal energy storage capacity [kWt$\cdot$h]
        param_dict['Lr']     = self.Lr.to('')          #L^r: Receiver pumping power per unit power produced [kWe/kWt]
        param_dict['Qrl']    = self.Qrl.to('kWh')      #Q^{rl}: Minimum operational thermal power delivered by receiver [kWt$\cdot$h]
        param_dict['Qrsb']   = self.Qrsb.to('kWh')     #Q^{rsb}: Required thermal power for receiver standby [kWt$\cdot$h]
        param_dict['Qrsd']   = self.Qrsd.to('kWh')     #Q^{rsd}: Required thermal power for receiver shut down [kWt$\cdot$h] 
        param_dict['Qru']    = self.Qru.to('kW')       #Q^{ru}: Allowable power per period for receiver start-up [kWt$\cdot$h]
        param_dict['Wh']     = self.Wh.to('kW')        #W^h: Heliostat field tracking parasitic loss [kWe]
        param_dict['Wnht']   = self.Wnht.to('kW')      #W^{ht}: Tower piping heat trace parasitic loss [kWe]
        
        return param_dict