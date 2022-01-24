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

class NuclearDispatch(GeneralDispatch):
    """
    The NuclearDispatch class is meant to set up and run Dispatch
    optimization as a mixed integer linear program problem using Pyomo,
    specifically for the NuclearTES NE2+SSC module.
    """
    
    def __init__(self, params, unitRegistry):
        """ Initializes the NuclearDispatch module
        
        The instantiation of this class receives a parameter dictionary from
        the NE2 module (created using the NuclearDispatchWrapper class). It calls
        on the GeneralDispatch __init__ to create the model. The GeneralDispatcher first
        creates an empty Concrete Model from Pyomo, then generates Parameters
        from the parameter dictionary, Variables, Objectives and Constraints.
        
        Inputs:
            params (dict)                : dictionary of Pyomo dispatch parameters
            unitRegistry (pint.registry) : unique unit Pint unit registry
        """
        
        # initialize Generic module, csv data arrays should be saved here
        GeneralDispatch.__init__( self, params, unitRegistry )


    def generate_params(self, params, skip_parent=False):
        """ Method to generate parameters within Pyomo Nuclear Model
        
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
        if not skip_parent:
            GeneralDispatch.generate_params(self, params)
        
        # lambdas to convert units and data to proper syntax
        gd = self.gd
        gu = self.gu 
        
        ### Cost Parameters ### 
        self.model.Cnuc = pe.Param(mutable=True, initialize=gd("Cnuc"), units=gu("Cnuc"))              #C^{nuc}: Operating cost of nuclear plant [\$/kWt$\cdot$h]
        self.model.Cnsu = pe.Param(mutable=True, initialize=gd("Cnsu"), units=gu("Cnsu"))              #C^{nsu}: Penalty for nuclear cold start-up [\$/start]
        self.model.Cnhsp = pe.Param(mutable=True, initialize=gd("Cnhsp"), units=gu("Cnhsp"))           #C^{nhsp}: Penalty for nuclear hot start-up [\$/start]
        
        ### Nuclear Parameters ###
        self.model.deltanl = pe.Param(mutable=True, initialize=gd("deltanl"), units=gu("deltanl")) #\delta^{nl}: Minimum time to start the nuclear plant [hr]
        self.model.En = pe.Param(mutable=True, initialize=gd("En"), units=gu("En"))                #E^n: Required energy expended to start nuclear plant [kWt$\cdot$h]
        self.model.Ln = pe.Param(mutable=True, initialize=gd("Ln"), units=gu("Ln"))                #L^n: Nuclear pumping power per unit power produced [kWe/kWt]
        self.model.Qnl = pe.Param(mutable=True, initialize=gd("Qnl"), units=gu("Qnl"))             #Q^{nl}: Minimum operational thermal power delivered by nuclear [kWt$\cdot$h]
        self.model.Qnsb = pe.Param(mutable=True, initialize=gd("Qnsb"), units=gu("Qnsb"))          #Q^{nsb}: Required thermal power for nuclear standby [kWt$\cdot$h]
        self.model.Qnsd = pe.Param(mutable=True, initialize=gd("Qnsd"), units=gu("Qnsd"))          #Q^{nsd}: Required thermal power for nuclear shut down [kWt$\cdot$h] 
        self.model.Qnu = pe.Param(mutable=True, initialize=gd("Qnu"), units=gu("Qnu"))             #Q^{nu}: Allowable power per period for nuclear start-up [kWt$\cdot$h]
        self.model.Wnht = pe.Param(mutable=True, initialize=gd("Wnht"), units=gu("Wnht"))          #W^{nht}: Nuclear piping heat trace parasitic loss [kWe]

        ### Time series Nuclear Parameters ###
        self.model.delta_ns = pe.Param(self.model.T, mutable=True, initialize=gd("delta_ns"), units=gu("delta_ns"))     #\delta^{ns}_{t}: Estimated fraction of period $t$ required for nuclear start-up [-]
        self.model.Qin_nuc = pe.Param(self.model.T, mutable=True, initialize=gd("Qin_nuc"), units=gu("Qin_nuc"))        #Q^{in,nuc}_{t}: Available thermal power generated by the nuclear plant in period $t$ [kWt]
        
        ### Initial Condition Parameters ###
        self.model.unsu0 = pe.Param(mutable=True, initialize=gd("unsu0"), units=gu("unsu0"))  #u^{nsu}_0: Initial nuclear start-up energy inventory [kWt$\cdot$h]
        self.model.yn0 = pe.Param(mutable=True, initialize=gd("yn0"), units=gu("yn0"))        #y^n_0: 1 if nuclear is generating ``usable'' thermal power initially, 0 otherwise  [az] this is new.
        self.model.ynsb0 = pe.Param(mutable=True, initialize=gd("ynsb0"), units=gu("ynsb0"))  #y^{nsb}_0: 1 if nuclear is in standby mode initially, 0 otherwise [az] this is new.
        self.model.ynsu0 = pe.Param(mutable=True, initialize=gd("ynsu0"), units=gu("ynsu0"))  #y^{nsu}_0: 1 if nuclear is in starting up initially, 0 otherwise [az] this is new.
        

    def generate_variables(self, skip_parent=False):
        """ Method to generate parameters within Pyomo Nuclear Model
        
        This method instantiates variables for the Pyomo Concrete Model, with
        domains. Does not need initial guesses here, they are defined in the 
        parameters. We first define continuous and binary variables for the 
        Power Cycle through GeneralDispatch, then declare nuclear variables.
        """
        
        # generating GeneralDispatch variables first (PowerCycle, etc.)
        if not skip_parent:
            GeneralDispatch.generate_variables(self)
        
        u = self.u_pyomo
        ### Decision Variables ###
        #------- Variables ---------
        self.model.unsu = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kWh)     #u^{nsu}: Nuclear start-up energy inventory at period $t$ [kWt$\cdot$h]
        self.model.xn = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)	    #x^n: Thermal power delivered by the nuclear at period $t$ [kWt]
        self.model.xnsu = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)      #x^{nsu}: Nuclear start-up power consumption at period $t$ [kWt]
        
        #------- Binary Variables ---------
        self.model.yn = pe.Var(self.model.T, domain=pe.Binary)        #y^n: 1 if nuclear is generating ``usable'' thermal power at period $t$; 0 otherwise
        self.model.ynhsp = pe.Var(self.model.T, domain=pe.Binary)	  #y^{nhsp}: 1 if nuclear hot start-up penalty is incurred at period $t$ (from standby); 0 otherwise
        self.model.ynsb = pe.Var(self.model.T, domain=pe.Binary)	  #y^{nsb}: 1 if nuclear is in standby mode at period $t$; 0 otherwise
        self.model.ynsd = pe.Var(self.model.T, domain=pe.Binary)	  #y^{nsd}: 1 if nuclear is shut down at period $t$; 0 otherwise
        self.model.ynsu = pe.Var(self.model.T, domain=pe.Binary)      #y^{nsu}: 1 if nuclear is starting up at period $t$; 0 otherwise
        self.model.ynsup = pe.Var(self.model.T, domain=pe.Binary)     #y^{nsup}: 1 if nuclear cold start-up penalty is incurred at period $t$ (from off); 0 otherwise


    def add_objective(self):
        """ Method to add an objective function to the Pyomo Nuclear Model
        
        This method adds an objective function to the Pyomo Nuclear Dispatch
        model. Typically, the LORE team defined a nested function and passes it
        into the Pyomo model. 
        """
        def objectiveRule(model):
            """ Maximize profits vs. costs """
            return (
                    sum( model.D[t] * 
                    #obj_profit
                    model.Delta[t]*model.P[t]*0.5*(model.wdot_s[t] - model.wdot_p[t])
                    #obj_cost_cycle_su_hs_sd
                    - (0.5*model.Ccsu*model.ycsup[t] + 0.2*model.Cchsp*model.ychsp[t] + model.alpha*model.ycsd[t])
                    #obj_cost_cycle_ramping
                    - model.Delta[t]*(model.C_delta_w*(model.wdot_delta_plus[t]+model.wdot_delta_minus[t])+model.C_v_w*(model.wdot_v_plus[t] + model.wdot_v_minus[t]))
                    #obj_cost_nuc_su_hs_sd
                    - (model.Cnsu*model.ynsup[t] + model.Cnhsp*model.ynhsp[t] + model.alpha*model.ynsd[t])
                    #obj_cost_ops
                    - model.Delta[t]*model.Cpc*model.wdot[t] - model.Delta[t]*model.Ccsb*model.Qb*model.ycsb[t] - model.Delta[t]*model.Cnuc*model.xn[t]
                    for t in model.T) 
                    )
        
        self.model.OBJ = pe.Objective(rule=objectiveRule, sense = pe.maximize)


    def addNuclearStartupConstraints(self):
        """ Method to add nuclear startup constraints to the Pyomo Nuclear Model
        
        This method adds constraints pertaining to nuclear reactor startup within the 
        Pyomo Nuclear Dispatch class. Several nested functions are defined. We don't
        model nuclear reactor startup in SSC, this is kept for completeness' sake. This
        is essentially used as a sanity check, if Pyomo deems that nuclear should shutdown
        or startup then something is wrong and SSC will complain. 
        """
        def nuc_inventory_rule(model,t):
            """ NP energy inventory for startup is balanced """
            if t == 1:
                return model.unsu[t] <= model.unsu0 + model.Delta[t]*model.xnsu[t]
            return model.unsu[t] <= model.unsu[t-1] + model.Delta[t]*model.xnsu[t]
        def nuc_inv_nonzero_rule(model,t):
            """ NP energy inventory is positive if starting up """
            return model.unsu[t] <= model.En * model.ynsu[t]
        def nuc_startup_rule(model,t):
            """ NP resumes normal operation after startup or if previously operating """
            if t == 1:
                return model.yn[t] <= model.unsu[t]/model.En + model.yn0 + model.ynsb0
            return model.yn[t] <= model.unsu[t]/model.En + model.yn[t-1] + model.ynsb[t-1]
        def nuc_su_persist_rule(model,t):
            """ NP cannot startup if operating in previous timestep """
            if t == 1: 
                return model.ynsu[t] + model.yn0 <= 1
            return model.ynsu[t] +  model.yn[t-1] <= 1
        def ramp_limit_rule(model,t):
            """ NP must abide by ramp rate limits """
            return model.xnsu[t] <= model.Qnu*model.ynsu[t]
        def nontrivial_nuclear_rule(model,t):
            """ NP trivial power prevents startup """
            return model.ynsu[t] <= model.Qin_nuc[t] / model.Qnl
        
        self.model.nuc_inventory_con = pe.Constraint(self.model.T,rule=nuc_inventory_rule)
        self.model.nuc_inv_nonzero_con = pe.Constraint(self.model.T,rule=nuc_inv_nonzero_rule)
        self.model.nuc_startup_con = pe.Constraint(self.model.T,rule=nuc_startup_rule)
        self.model.nuc_su_persist_con = pe.Constraint(self.model.T,rule=nuc_su_persist_rule)
        self.model.nuc_ramp_limit_con = pe.Constraint(self.model.T,rule=ramp_limit_rule)
        self.model.nontrivial_nuclear_con = pe.Constraint(self.model.T,rule=nontrivial_nuclear_rule)
        
        
    def addNuclearSupplyAndDemandConstraints(self):
        """ Method to add nuclear supply and demand constraints to the Pyomo Nuclear Model
        
        This method adds constraints pertaining to nuclear supply and demand energy
        constraints. Some constraints might be redundant, they are adapted from the CSP
        constraints (thanks LORE team).
        """
        def nuc_production_rule(model,t):
            """ Upper bound on thermal energy produced by NP """
            return model.xn[t] + model.xnsu[t] + model.Qnsd*model.ynsd[t] <= model.Qin_nuc[t]
        def nuc_generation_rule(model,t):
            """ Thermal energy production by NP only when operating """
            return model.xn[t] <= model.Qin_nuc[t] * model.yn[t]
        def min_generation_rule(model,t):
            """ Lower bound on thermal energy produced by NP """
            return model.xn[t] >= model.Qnl * model.yn[t]
        def nuc_gen_persist_rule(model,t):
            """ NP not able to operate if no thermal power (adapted from CSP) """
            return model.yn[t] <= model.Qin_nuc[t] / model.Qnl 
        
        self.model.nuc_production_con = pe.Constraint(self.model.T,rule=nuc_production_rule)
        self.model.nuc_generation_con = pe.Constraint(self.model.T,rule=nuc_generation_rule)
        self.model.nuc_min_generation_con = pe.Constraint(self.model.T,rule=min_generation_rule)
        self.model.nuc_gen_persist_con = pe.Constraint(self.model.T,rule=nuc_gen_persist_rule)
        
        
    def addNuclearNodeLogicConstraints(self):
        """ Method to add nuclear logic constraints to the Pyomo General Model
        
        This method adds constraints pertaining to tracking binary variable
        logic for the nuclear plant. Essentially, to make sure the correct modes
        are activated when they are allowed. There is no Nuclear standby (nor CSP)
        so those particular constraints are also just for completeness sake. 
        """
        def nuc_su_sb_persist_rule(model,t):
            """ NP startup and standby cannot coincide """
            return model.ynsu[t] + model.ynsb[t] <= 1
        def nuc_sb_persist_rule(model,t):
            """ NP standby and operation cannot coincide """
            return model.yn[t] + model.ynsb[t] <= 1
        def rsb_persist_rule(model,t):
            """ NP standby must follow operating or another standby timestep """
            if t == 1:
                return model.ynsb[t] <= (model.yn0 + model.ynsb0) 
            return model.ynsb[t] <= model.yn[t-1] + model.ynsb[t-1]
        def nuc_su_pen_rule(model,t):
            """ NP cold startup penalty """
            if t == 1:
                return model.ynsup[t] >= model.ynsu[t] - model.ynsu0 
            return model.ynsup[t] >= model.ynsu[t] - model.ynsu[t-1]
        def nuc_hs_pen_rule(model,t):
            """ NP hot startup penalty: if we went from standby to operating """
            if t == 1:
                return model.ynhsp[t] >= model.yn[t] - (1 - model.ynsb0)
            return model.ynhsp[t] >= model.yn[t] - (1 - model.ynsb[t-1])
        def nuc_shutdown_rule(model,t):
            """ NP shutdown protocol """
            current_Delta = model.Delta[t]
            # structure of inequality is lb <= model.param <= ub with strict=False by default
            if self.eval_ineq(1,current_Delta) and t == 1: #not strict
                return 0 >= model.yn0 - model.yn[t] +  model.ynsb0 - model.ynsb[t]
            elif self.eval_ineq(1,current_Delta) and t > 1: # not strict
                return model.ynsd[t-1] >= model.yn[t-1] - model.yn[t] + model.ynsb[t-1] - model.ynsb[t]
            elif self.eval_ineq(current_Delta,1,strict=True) and t == 1:
                return model.ynsd[t] >= model.yn0  - model.yn[t] + model.ynsb0 - model.ynsb[t]
            # only case remaining: Delta[t]<1, t>1
            return model.ynsd[t] >= model.yn[t-1] - model.yn[t] + model.ynsb[t-1] - model.ynsb[t]
        
        self.model.nuc_su_sb_persist_con = pe.Constraint(self.model.T,rule=nuc_su_sb_persist_rule)
        self.model.nuc_sb_persist_con = pe.Constraint(self.model.T,rule=nuc_sb_persist_rule)
        self.model.nuc_rsb_persist_con = pe.Constraint(self.model.T,rule=rsb_persist_rule)
        self.model.nuc_su_pen_con = pe.Constraint(self.model.T,rule=nuc_su_pen_rule)
        self.model.nuc_hs_pen_con = pe.Constraint(self.model.T,rule=nuc_hs_pen_rule)
        self.model.nuc_shutdown_con = pe.Constraint(self.model.T,rule=nuc_shutdown_rule)
        

    def addTESEnergyBalanceConstraints(self):
        """ Method to add TES constraints to the Pyomo Nuclear Model
        
        This method adds constraints pertaining to TES energy balance from charging
        with thermal power and discharging to the power cycle. 
        
        TODO: revisit tes_start_up_rule -> do we need this for nuclear?
        TODO: do we need maintain_tes_rule?
        """
        def tes_balance_rule(model, t):
            """ Balance of energy to and from TES """
            if t == 1:
                return model.s[t] - model.s0 == model.Delta[t] * (model.xn[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.x[t] + model.Qnsb*model.ynsb[t]))
            return model.s[t] - model.s[t-1] == model.Delta[t] * (model.xn[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.x[t] + model.Qnsb*model.ynsb[t]))
        def tes_upper_rule(model, t):
            """ Upper bound to TES charge state """
            return model.s[t] <= model.Eu
        def tes_start_up_rule(model, t):
            """ Ensuring sufficient TES charge level to startup NP """
            if t == 1:
                return model.s0 >= model.Delta[t]*model.delta_ns[t]*( (model.Qu + model.Qb)*( -3 + model.ynsu[t] + model.y0 + model.y[t] + model.ycsb0 + model.ycsb[t] ) + model.x[t] + model.Qb*model.ycsb[t] )
            return model.s[t-1] >= model.Delta[t]*model.delta_ns[t]*( (model.Qu + model.Qb)*( -3 + model.ynsu[t] + model.y[t-1] + model.y[t] + model.ycsb[t-1] + model.ycsb[t] ) + model.x[t] + model.Qb*model.ycsb[t] )
        def maintain_tes_rule(model):
            """ Final state of TES has to be less than or equal to start """
            return model.s[model.num_periods] <= model.s0
        
        self.model.tes_balance_con = pe.Constraint(self.model.T,rule=tes_balance_rule)
        self.model.tes_upper_con = pe.Constraint(self.model.T,rule=tes_upper_rule)
        self.model.tes_start_up_con = pe.Constraint(self.model.T,rule=tes_start_up_rule)
        self.model.maintain_tes_con = pe.Constraint(rule=maintain_tes_rule)


    def addPiecewiseLinearEfficiencyConstraints(self):
        """ Method to add efficiency constraints to the Pyomo Nuclear Model
        
        This method adds constraints pertaining to efficiency constraints defined
        as a piecewise linear approximation. Also referred to as Cycle supply and 
        demand constraints. In the NuclearDispatch, we add an extra balance of power
        with respect to energy storage and power produced from the nuclear reactor. 
        
        TODO: This should be revisited when adding MED!!
        """
        def grid_sun_rule(model, t):
            """ Balance of power flow, i.e. sold vs purchased """
            return (
                    model.wdot_s[t] - model.wdot_p[t] == (1-model.etac[t])*model.wdot[t]
                		- model.Ln*(model.xn[t] + model.xnsu[t] + model.Qnl*model.ynsb[t])
                		- model.Lc*model.x[t] 
                        - model.Wb*model.ycsb[t] - model.Wnht*(model.ynsb[t]+model.ynsu[t])		#Is Wrsb energy [kWh] or power [kW]?  [az] Wrsb = Wht in the math?
            )
        
        # call the parent version of this method
        GeneralDispatch.addPiecewiseLinearEfficiencyConstraints(self)
        
        # additional constraints
        self.model.grid_sun_con = pe.Constraint(self.model.T,rule=grid_sun_rule)


    def generate_constraints(self, skip_parent=False):
        """ Method to add ALL constraints to the Pyomo Nuclear Model
        
        This method calls the previously defined constraint methods to instantiate
        them and add to the existing model. This method first calls the GeneralDispatch
        version to set PowerCycle constraints, then calls nuclear constraint methods
        to add them to the model. 
        """
        
        # generating GeneralDispatch constraints first (PowerCycle, etc.)
        if not skip_parent:
            # call general PC constraints
            GeneralDispatch.generate_constraints(self)
            # call TES energy balance constraints specific to current dispatch model
            self.addTESEnergyBalanceConstraints()
        
        self.addNuclearStartupConstraints()
        self.addNuclearSupplyAndDemandConstraints()
        self.addNuclearNodeLogicConstraints()
        

    
# =============================================================================
# Dispatch Wrapper
# =============================================================================

class NuclearDispatchParamWrap(GeneralDispatchParamWrap):
    """
    The NuclearDispatchParamWrap class is meant to be the staging area for the 
    creation of Parameters ONLY for the NuclearDispatch class. It communicates 
    with the NE2 modules, receiving SSC and PySAM input dictionaries to calculate 
    both static parameters used for every simulation segment AND initial conditions 
    that can be updated.
    """
    
    def __init__(self, unit_registry, SSC_dict=None, PySAM_dict=None, pyomo_horizon=48, 
                   dispatch_time_step=1):
        """ Initializes the NuclearDispatchParamWrap module
        
        Inputs:
            unitRegistry (pint.registry)   : unique unit Pint unit registry
            SSC_dict (dict)                : dictionary of SSC inputs needed to run modules
            PySAM_dict (dict)              : dictionary of PySAM inputs + file names
            pyomo_horizon (int Quant)      : length of Pyomo simulation segment (hours)
            dispatch_time_step (int Quant) : length of each Pyomo time step (hours)
        """
        
        GeneralDispatchParamWrap.__init__( self, unit_registry, SSC_dict, PySAM_dict, 
                            pyomo_horizon, dispatch_time_step )


    def set_design(self, skip_parent=False):
        """ Method to calculate and save design point values of Plant operation
        
        This method extracts values and calculates for design point parameters 
        of our Plant (e.g., nuclear thermal power output, power cycle efficiency,
        inlet and outlet temperatures, etc.). 
        """
        
        u = self.u
        
        if not skip_parent:
            GeneralDispatchParamWrap.set_design(self)
        
        # nuclear parameters
        self.q_nuc_design = self.SSC_dict['q_dot_nuclear_des'] * u.MW      # nuclear design thermal power
        
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
        
        
    def set_fixed_cost_parameters(self, param_dict, skip_parent=False):
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
        if not skip_parent:
            param_dict = GeneralDispatchParamWrap.set_fixed_cost_parameters( self, param_dict )
        
        C_nuc  = self.PySAM_dict['nuc_op_cost'] * u.USD / u.MWh  # value taken from Cory @ Westinghouse, later converted to $/kWh        
        C_nsu  = self.PySAM_dict['nuc_cold_su'] * u.USD
        C_nhsp = self.PySAM_dict['nuc_hot_su'] * u.USD
        
        # original nuc_op_cost given in MWh-electric rather than MWh-thermal
        C_nuc *= self.eta_design # eta => electric/thermal conversion
        
        ### Cost Parameters ###
        param_dict['Cnuc']   = C_nuc.to('USD/kWh')  #C^{nuc}: Operating cost of nuclear plant [\$/kWt$\cdot$h]
        param_dict['Cnsu']   = C_nsu.to('USD')      #C^{nsu}: Penalty for nuclear cold start-up [\$/start]
        param_dict['Cnhsp']  = C_nhsp.to('USD')     #C^{nhsp}: Penalty for nuclear hot start-up [\$/start]

        return param_dict
        

    def set_nuclear_parameters(self, param_dict):
        """ Method to set parameters specific to the Nuclear Plant for Dispatch optimization
        
        This method calculates some parameters specific to the NuclearTES plant
        which are meant to be fixed throughout the simulation. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        # grabbing unit registry set up in GeneralDispatch
        u = self.u 
        
        dw_rec_pump             = self.PySAM_dict['dw_nuc_pump']*u.MW   # TODO: Pumping parasitic at design point reciever mass flow rate (MWe)
        tower_piping_ht_loss    = self.PySAM_dict['nuc_piping_ht_loss']*u.kW   # TODO: Tower piping heat trace full-load parasitic load (kWe) 
        q_nuc_standby_fraction  = self.PySAM_dict['q_nuc_standby_frac']        # TODO: Nuclear standby energy consumption (fraction of design point thermal power)
        q_nuc_shutdown_fraction = self.PySAM_dict['q_nuc_shutdown_frac']       # TODO: Nuclear shutdown energy consumption (fraction of design point thermal power)
        
        self.deltanl = self.SSC_dict['rec_su_delay']*u.hr
        self.En     = (self.SSC_dict['rec_qf_delay'] * u.kWh / u.kW ) * self.q_nuc_design 
        self.Eu     = self.SSC_dict['tshours']*u.hr * self.q_pb_design
        self.Ln     = dw_rec_pump / self.q_nuc_design
        self.Qnl    = self.SSC_dict['f_nuc_min'] * self.q_nuc_design 
        self.Qnsb   = q_nuc_standby_fraction  * self.q_nuc_design 
        self.Qnsd   = q_nuc_shutdown_fraction * self.q_nuc_design 
        self.Qnu    = self.En / self.deltanl  
        self.Wnht    = tower_piping_ht_loss
        
        ### CSP Field and Receiver Parameters ###
        param_dict['deltanl'] = self.deltanl.to('hr')  #\delta^{nl}: Minimum time to start the nuclear plant [hr]
        param_dict['En']     = self.En.to('kWh')       #E^n: Required energy expended to start nuclear plant [kWt$\cdot$h]
        param_dict['Eu']     = self.Eu.to('kWh')       #E^u: Thermal energy storage capacity [kWt$\cdot$h]
        param_dict['Ln']     = self.Ln.to('')          #L^n: Nuclear pumping power per unit power produced [kWe/kWt]
        param_dict['Qnl']    = self.Qnl.to('kW')      #Q^{nl}: Minimum operational thermal power delivered by nuclear [kWt$\cdot$h]
        param_dict['Qnsb']   = self.Qnsb.to('kW')     #Q^{nsb}: Required thermal power for nuclear standby [kWt$\cdot$h]
        param_dict['Qnsd']   = self.Qnsd.to('kW')     #Q^{nsd}: Required thermal power for nuclear shut down [kWt$\cdot$h] 
        param_dict['Qnu']    = self.Qnu.to('kW')       #Q^{nu}: Allowable power per period for nuclear start-up [kWt$\cdot$h]
        param_dict['Wnht']   = self.Wnht.to('kW')      #W^{nht}: Nuclear piping heat trace parasitic loss [kWe]
        
        return param_dict
    
    
    def set_time_series_nuclear_parameters(self, param_dict, updated_dict=None):
        """ Method to set fixed costs of the Plant for Dispatch optimization
        
        This method calculates some time series parameters for the Plant operations, startup,
        standby, etc. These are NOT meant to be fixed, but updated at the beginning
        of every segment using the latest SSC outputs or to extract the next relevant
        segment of pricing arrays, efficiencies, etc. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
            updated_dict (dict) : dictionary with updated SSC initial conditions from previous run  
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        #MAKE SURE TO CALL THIS METHOD AFTER THE NUCLEAR PARAMETERS 
        u = self.u
        
        # differentiating between first and updated runs
        if updated_dict is None:
            self.current_Plant = copy.deepcopy(self.SSC_dict)
            self.first_run = True
        else:
            self.current_Plant = updated_dict
            self.first_run = False
        
        # thermal power and startup params
        self.Dnsu    = self.PySAM_dict['Dnsu']*u.hr   # Minimum time to start the nuclear plant (hr)
        self.Qin_nuc = self.current_Plant['Q_nuc_thermal'] * u.MW    # value here taken from a previous Plant-SSC run

        if len(self.Qin_nuc) < len(param_dict['Delta']):
            self.Qin_nuc = np.hstack( [self.current_Plant['Q_nuc_thermal'], self.current_Plant['Q_nuc_thermal']] ) * u.MW
        
        # instantiating arrays
        n  = len(self.Delta)
        delta_ns = np.zeros(n)
        
        # loop to set startup nuclear params (probably won't need, keeping it for unit testing?)
        for t in range(n):
            Ein = self.Qin_nuc[t]*self.Delta[t]
            E_compare = (self.En / max(1.*u.kWh, Ein.to('kWh'))).to('')
            delta_ns[t] = min(1., max( E_compare, self.Dnsu/self.Delta[t]))
        
        self.delta_ns   = delta_ns

        ### Time series CSP Parameters ###
        param_dict['delta_ns']  = self.delta_ns      #\delta^{ns}_{t}: Estimated fraction of period $t$ required for nuclear start-up [-]
        param_dict['Qin_nuc']   = self.Qin_nuc.to('kW')  #Q^{in}_{t}: Available thermal power generated by the nuclear plant in period $t$ [kWt]
        
        return param_dict


    def set_initial_state(self, param_dict, updated_dict=None, plant=None, npts=None, skip_parent=False):
        """ Method to set the initial state of the Plant before Dispatch optimization
        
        This method uses SSC data to set the initial state of the Plant before Dispatch
        optimization in Pyomo. This method is called in two ways: once before starting 
        the simulation loop, in which case it only uses values from the SSC_dict portion
        of the given JSON script. The method is also called within the simulation loop
        to update the initial state parameters based on the ending conditions of the 
        previous simulation segment (provided by SSC). 
        
        TODO: can we just input another dictionary instead of passing the full Plant?
        
        Inputs:
            param_dict (dict)    : dictionary of Pyomo dispatch parameters
            updated_dict (dict)  : dictionary with updated SSC initial conditions from previous run
            plant (obj)          : the full PySAM Plant object. 
            npts (int)           : length of the SSC horizon
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        u = self.u
        
        # First filling out initial states from GeneralDispatcher
        if not skip_parent:
            param_dict = GeneralDispatchParamWrap.set_initial_state( self, param_dict, updated_dict, plant, npts )
        
        if updated_dict is None:
            self.current_Plant = copy.deepcopy(self.SSC_dict)
            self.first_run = True
        else:
            self.current_Plant = updated_dict
            self.first_run = False
        
        # TES masses, temperatures, specific heat
        m_hot  = self.m_tes_design * (self.current_Plant['csp.pt.tes.init_hot_htf_percent']/100)        # Available active mass in hot tank
        T_tes_hot_init  = (self.current_Plant['T_tank_hot_init']*u.celsius).to('degK')
        T_tes_init  = 0.5*(T_tes_hot_init + self.T_htf_cold)
        cp_tes_init = SSCHelperMethods.get_cp_htf(self.u, T_tes_init, self.SSC_dict['rec_htf'] )
        
        # important parameters
        s_current          = m_hot * cp_tes_init * (T_tes_hot_init - self.T_htf_cold) # TES capacity
        s0                 = min(self.Eu.to('kWh'), s_current.to('kWh')  )
        yn0                = (self.current_Plant['nuc_op_mode_initial'] == 2)
        ynsb0              = False   # We don't have standby mode for either Nuclear or CSP
        ynsu0              = (self.current_Plant['nuc_op_mode_initial'] == 1)
        t_nuc              = self.current_Plant['nuc_startup_time_remain_init']
        t_nuc_suinitremain = t_nuc if not np.isnan( t_nuc ) else 0.0
        e_nuc              = self.current_Plant['nuc_startup_energy_remain_init']
        e_nuc_suinitremain = e_nuc if not np.isnan( e_nuc ) else 0.0
        nuc_accum_time     = max(0.0*u.hr, self.Dnsu - t_nuc_suinitremain*u.hr )
        nuc_accum_energy   = max(0.0*u.Wh, self.En   - e_nuc_suinitremain*u.Wh )

        # defining parameters
        self.s0    = s0              #s_0: Initial TES reserve quantity  [kWt$\cdot$h]
        self.yn0   = yn0             #y^r_0: 1 if nuclear is generating ``usable'' thermal power initially, 0 otherwise 
        self.ynsb0 = ynsb0           #y^{nsb}_0: 1 if nuclear is in standby mode initially, 0 otherwise
        self.ynsu0 = ynsu0           #y^{nsu}_0: 1 if nuclear is in starting up initially, 0 otherwise
        
        # Initial nuclear startup energy inventory
        self.unsu0 = min(nuc_accum_energy, nuc_accum_time * self.Qnu)  # Note, SS receiver model in ssc assumes full available power is used for startup (even if, time requirement is binding)
        if self.unsu0 > (1.0 - 1.e-6)*self.En:
            self.unsu0 = self.En

        param_dict['s0']     = self.s0.to('kWh')    #s_0: Initial TES reserve quantity  [kWt$\cdot$h]
        param_dict['unsu0']  = self.unsu0.to('kWh') #u^{nsu}_0: Initial nuclear start-up energy inventory [kWt$\cdot$h]
        param_dict['yn0']    = self.yn0             #y^n_0: 1 if nuclear is generating ``usable'' thermal power initially, 0 otherwise
        param_dict['ynsb0']  = self.ynsb0           #y^{nsb}_0: 1 if nuclear is in standby mode initially, 0 otherwise
        param_dict['ynsu0']  = self.ynsu0           #y^{nsu}_0: 1 if nuclear is in starting up initially, 0 otherwise
        
        return param_dict
    
    
# =============================================================================
# Dispatch Outputs
# =============================================================================
  
class NuclearDispatchOutputs(object):
    """
    The NuclearDispatchOutputs class is meant to handle outputs from a given,
    solved Pyomo Dispatch model. It returns desired outputs in appropriate formats
    and syntaxes for PostProcessing and linking simulation segments between Pyomo
    and SSC calls. 
    """
    
    def get_dispatch_targets_from_Pyomo(dispatch_model, horizon, N_full, run_loop=False):
        """ Method to set fixed costs of the Plant
        
        This method parses through the solved Pyomo model for Dispatch optimization
        and extracts results that are used as Dispatch Targets in the *SAME* simulation
        segment but in SSC rather than Pyomo. If we're not running a loop, we can
        still update SSC only I guess this happens once for whatever Pyomo horizon
        is defined (this might not be a feature we keep long-term, perhaps only for
                    debugging). 
        
        Inputs:
            dispatch_model (Pyomo model) : solved Pyomo Dispatch model (ConcreteModel)
            horizon (float Quant)        : length of time of horizon, whether SSC or Pyomo (in hours)
            N_full (int)                 : length of full simulation time (in hours, no Quant)
            run_loop (bool)              : flag to determine if simulation is segmented
        Outputs:
            disp_targs (dict) : dictionary of dispatch target arrays for use in SSC 
        """
        
        dm = dispatch_model
        
        # range of pyomo and SSC horizon times
        t_pyomo = dm.model.T
        f_ind   = int( horizon.to('hr').m ) # index in hours of final horizon (e.g. 24)
        t_horizon = range(f_ind)
        
        # if we're not running a loop, define a list of 0s to pad the output so it matches full array size
        if not run_loop:
            N_leftover = N_full - f_ind
            empty_array = [0]*N_leftover
        
        #----Receiver Binary Outputs----
        yn   = np.array([pe.value(dm.model.yn[t])   for t in t_pyomo])
        ynsu = np.array([pe.value(dm.model.ynsu[t]) for t in t_pyomo])
        ynsb = np.array([pe.value(dm.model.ynsb[t]) for t in t_pyomo])
        
        #----Cycle Binary Outputs----
        y    = np.array([pe.value(dm.model.y[t])    for t in t_pyomo])
        ycsu = np.array([pe.value(dm.model.ycsu[t]) for t in t_pyomo])
        ycsb = np.array([pe.value(dm.model.ycsb[t]) for t in t_pyomo])
    
        #----Cycle Thermal Power Utilization----
        x = np.array([pe.value(dm.model.x[t])   for t in t_pyomo])/1000. # from kWt -> MWt
        
        #----Thermal Capacity for Cycle Startup and Operation----
        Qc = np.array([pe.value(dm.model.Qc[t]) for t in t_pyomo])/1000. # from kWt -> MWt
        Qu = dm.model.Qu.value/1000. # from kWt -> MWt
    
        # dispatch target -- nuclear startup/standby binaries
        is_nuc_su_allowed_in = [1 if (yn[t] + ynsu[t] + ynsb[t]) > 0.001 else 0 for t in t_horizon]  # Nuclear on, startup, or standby
        is_nuc_sb_allowed_in = [1 if ynsb[t] > 0.001                     else 0 for t in t_horizon]  # Nuclear standby
        
        # dispatch target -- cycle startup/standby binaries
        is_pc_su_allowed_in  = [1 if (y[t] + ycsu[t]) > 0.001 else 0 for t in t_horizon]  # Cycle on or startup
        is_pc_sb_allowed_in  = [1 if ycsb[t] > 0.001          else 0 for t in t_horizon]  # Cycle standby
    
        # dispatch target -- cycle thermal inputs and capacities
        q_pc_target_su_in    = [Qc[t] if ycsu[t] > 0.001 else 0.0 for t in t_horizon]
        q_pc_target_on_in    = [x[t]                              for t in t_horizon]
        q_pc_max_in          = [Qu                                for t in t_horizon]
        
        # empty dictionary for output
        disp_targs = {}
        
        # if we're running full simulation in steps, save SSC horizon portion of Pyomo horizon results
        if run_loop:
            disp_targs['is_rec_su_allowed_in'] = is_nuc_su_allowed_in 
            disp_targs['is_rec_sb_allowed_in'] = is_nuc_sb_allowed_in
            disp_targs['is_pc_su_allowed_in']  = is_pc_su_allowed_in 
            disp_targs['is_pc_sb_allowed_in']  = is_pc_sb_allowed_in  
            disp_targs['q_pc_target_su_in']    = q_pc_target_su_in  
            disp_targs['q_pc_target_on_in']    = q_pc_target_on_in
            disp_targs['q_pc_max_in']          = q_pc_max_in
        # if we're running full simulation all at once, need arrays to match size of full sim
        # TODO: is this a feature we want in the long term? Or just for debugging the first Pyomo call?
        else:
            disp_targs['is_rec_su_allowed_in'] = np.hstack( [is_nuc_su_allowed_in , empty_array] ).tolist()
            disp_targs['is_rec_sb_allowed_in'] = np.hstack( [is_nuc_sb_allowed_in , empty_array] ).tolist()
            disp_targs['is_pc_su_allowed_in']  = np.hstack( [is_pc_su_allowed_in  , empty_array] ).tolist()
            disp_targs['is_pc_sb_allowed_in']  = np.hstack( [is_pc_sb_allowed_in  , empty_array] ).tolist()
            disp_targs['q_pc_target_su_in']    = np.hstack( [q_pc_target_su_in    , empty_array] ).tolist()
            disp_targs['q_pc_target_on_in']    = np.hstack( [q_pc_target_on_in    , empty_array] ).tolist()
            disp_targs['q_pc_max_in']          = np.hstack( [q_pc_max_in          , empty_array] ).tolist()
            
        return disp_targs
    

