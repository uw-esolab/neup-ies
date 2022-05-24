#!/usr/bin/env python3
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
from dispatch.NuclearDispatch import NuclearDispatchParamWrap
from dispatch.NuclearDispatch import NuclearDispatch
import numpy as np
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
import os, copy

class IndirectNuclearDispatch(NuclearDispatch):
    """
    The IndirectNuclearDispatch class is meant to set up and run Dispatch
    optimization as a mixed integer linear program problem using Pyomo,
    specifically for the IndirectNuclearTES NE2+SSC module.
    """
    
    def __init__(self, **kwargs):
        """ Initializes the IndirectNuclearDispatch module
        
        The instantiation of this class receives a parameter dictionary from
        the NE2 module (created using the NuclearDispatchWrapper class). It calls
        on the GeneralDispatch __init__ to create the model. The GeneralDispatcher first
        creates an empty Concrete Model from Pyomo, then generates Parameters
        from the parameter dictionary, Variables, Objectives and Constraints.
        
        Inputs:
            params (dict)                : dictionary of Pyomo dispatch parameters
            unitRegistry (pint.registry) : unique unit Pint unit registry
        """
        
        # overriding base class default value unless we get a keyword from higher up
        kwargs['direct'] = kwargs['direct'] if 'direct' in kwargs else False
        super().__init__(**kwargs)


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
        
        super().generate_params(params)

        # lambdas to convert units and data to proper syntax
        gd = self.gd
        gu = self.gu 
        
        # New Cycle Parameters
        self.model.Wnc    = pe.Param(mutable=True, initialize=gd("Wnc"), units=gu("Wnc"))            #W^nc: Cycle capacity for nuclear power only [kWe]
        self.model.Qnhx   = pe.Param(mutable=True, initialize=gd("Qnhx"), units=gu("Qnhx"))       #Q^{nhx}: Upper limit on nuclear power to charge TES due to HX efficiency [kWt]
        self.model.eta_LD = pe.Param(mutable=True, initialize=gd("eta_LD"), units=gu("eta_LD"))   #\eta^{LD}: Linearized cycle efficiency during low demand operation [-]
        self.model.eta_HD = pe.Param(mutable=True, initialize=gd("eta_HD"), units=gu("eta_HD"))   #\eta^{HD}: Linearized cycle efficiency during high demand operation [-]
        self.model.intcpts_LD = pe.Param(mutable=True, initialize=gd("intcpts_LD"), units=gu("intcpts_LD"))   #\c^{LD}: Linearized cycle efficiency y-intercepts during low demand operation [-]
        self.model.intcpts_HD = pe.Param(mutable=True, initialize=gd("intcpts_HD"), units=gu("intcpts_HD"))   #\c^{HD}: Linearized cycle efficiency y-intercepts during high demand operation [-]
                
        
    def generate_variables(self, skip_parent=False):
        """ Method to generate parameters within Pyomo Nuclear Model
        
        This method instantiates variables for the Pyomo Concrete Model, with
        domains. Does not need initial guesses here, they are defined in the 
        parameters. We first define continuous and binary variables for the 
        Power Cycle through GeneralDispatch, then declare nuclear variables.
        """
        
        super().generate_variables()

        u = self.u_pyomo
        
        ### Decision Variables ###
        #------- Variables ---------
        self.model.xnp = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)    #\dot{x}^{np}_t: Thermal power delivered to cycle from nuclear at period $t$ [kWt]                            #x: Cycle thermal power utilization at period $t$ [kWt]
        self.model.xtesp = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)  #\dot{x}^{tesp}_t: Thermal power delivered to cycle from storage at period $t$ [kWt]
        self.model.xntes = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)  #\dot{x}^{ntes}_t: Thermal power delivered to storage from nuclear at period $t$ [kWt]

        #------- Binary Variables ---------
        self.model.ytesp = pe.Var(self.model.T, domain=pe.Binary)         #y: 1 if cycle is receiving thermal power from TES at period $t$; 0 otherwise
        self.model.yntes = pe.Var(self.model.T, domain=pe.Binary)         #y: 1 if storage is receiving thermal power from nuclear at period $t$; 0 otherwise


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
                    model.Delta[t]*model.P[t]*(model.wdot_s[t] - model.wdot_p[t])
                    #obj_cost_cycle_su_hs_sd
                    - (model.Ccsu*model.ycsup[t] + model.Cchsp*model.ychsp[t] + model.alpha*model.ycsd[t])
                    #obj_cost_cycle_ramping
                    - model.Delta[t]*(model.C_delta_w*(model.wdot_delta_plus[t]+model.wdot_delta_minus[t]) + model.C_v_w*(model.wdot_v_plus[t] + model.wdot_v_minus[t]))
                    #obj_cost_nuc_su_hs_sd
                    - (model.Cnsu*model.ynsup[t] + model.Cnhsp*model.ynhsp[t] + model.alpha*model.ynsd[t])
                    #obj_cost_ops
                    - model.Delta[t]*model.Cpc*model.wdot[t] - model.Delta[t]*model.Ccsb*model.Qb*model.ycsb[t] - model.Delta[t]*model.Cnuc*(model.xnp[t] + model.xntes[t])
                    for t in model.T) 
                    )
        
        self.model.OBJ = pe.Objective(rule=objectiveRule, sense = pe.maximize)
        
        
    def addPiecewiseLinearEfficiencyConstraints(self):
        """ Method to add efficiency constraints to the Pyomo IndirectNuclear Model
        
        This method adds constraints pertaining to efficiency constraints defined
        as a piecewise linear approximation. Also referred to as Cycle supply and 
        demand constraints. In the IndirectNuclearDispatch, we add an extra balance of power
        with respect to energy storage and power producfrom scipy.interpolate import interp1ded from the nuclear reactor. 
        
        TODO: This should be revisited when adding MED!!
        """
        
        super(NuclearDispatch, self).addPiecewiseLinearEfficiencyConstraints()
        
        def power_rule(model, t):
            """ Model of electric power vs. heat input as linear function """
            return model.wdot[t] <= (model.etaamb[t]/model.eta_des)*(model.eta_LD*(model.xnp[t] + model.xtesp[t]) 
                                                                     + model.y[t]*model.intcpts_LD)
        def high_demand_power_rule(model, t):
            """ NEW: Model of electric power vs. heat input as linear function for high demand efficiency """
            return model.wdot[t] <= (model.etaamb[t]/model.eta_des)*(model.eta_HD*(model.xnp[t] + model.xtesp[t])
                                                                    + model.y[t]*model.intcpts_HD)
        def grid_sun_rule(model, t):
            """ Balance of power flow, i.e. sold vs purchased """
            return (model.wdot_s[t] - model.wdot_p[t] == (1-model.etac[t])*model.wdot[t]
                		- model.Ln*(model.xnp[t] + model.xntes[t] + model.xnsu[t] + model.Qnl*model.ynsb[t])
                		- model.Lc*(model.xnp[t] + model.xtesp[t])
                        - model.Wb*model.ycsb[t] - model.Wnht*(model.ynsb[t]+model.ynsu[t])   )		#Is Wrsb energy [kWh] or power [kW]?  [az] Wrsb = Wht in the math?
        
        # overloaded constraints
        self.model.del_component( self.model.power_con )
        self.model.power_con = pe.Constraint(self.model.T,rule=power_rule)
        self.model.grid_sun_con = pe.Constraint(self.model.T,rule=grid_sun_rule)
        
        # new constraint
        self.model.hd_power_con = pe.Constraint(self.model.T,rule=high_demand_power_rule)


    def addCycleStartupConstraints(self):
        """ Method to add cycle startup constraints to the Pyomo General Model
        
        This method adds constraints pertaining to cycle startup within the Pyomo
        General Dispatch class. Several nested functions are defined. 
        """
        
        super(NuclearDispatch, self).addCycleStartupConstraints()
        
        def pc_production_rule(model, t):
            """ Heat input used for PC startup doesn't exceed max """
            return model.xnp[t] + model.xtesp[t] + model.Qc[t]*model.ycsu[t] <= model.Qu
        def pc_generation_rule(model, t):
            """ Upper bound on heat input to PC """
            return model.xnp[t] + model.xtesp[t] <= model.Qu * model.y[t]
        def pc_min_gen_rule(model, t):
            """ Lower bound on heat input to PC """
            return model.xnp[t] + model.xtesp[t] >= model.Ql * model.y[t]
        def pc_generation_thru_tes_rule(model, t):
            """ NEW: Upper bound on heat input from TES to PC """
            return model.xtesp[t] <= model.Qu * model.ytesp[t]
        def tes_dispatch_to_pc_rule(model, t):
            """ NEW: Enabling TES dispatch only when demand is higher than nominal """
            return model.ytesp[t] <= model.wdot[t] / model.Wnc
        def tes_dispatch_persist_rule(model, t):
            """ NEW: TES Dispatch only allowed if PC is on """
            return model.ytesp[t] <= model.y[t]
        
        # overloaded constraints
        self.model.del_component( self.model.pc_production_con )
        self.model.del_component( self.model.pc_generation_con )
        self.model.del_component( self.model.pc_min_gen_con )
        
        self.model.pc_production_con = pe.Constraint(self.model.T,rule=pc_production_rule)
        self.model.pc_generation_con = pe.Constraint(self.model.T,rule=pc_generation_rule)
        self.model.pc_min_gen_con = pe.Constraint(self.model.T,rule=pc_min_gen_rule)
        
        # new constraints
        self.model.pc_generation_thru_tes_con = pe.Constraint(self.model.T,rule=pc_generation_thru_tes_rule)
        self.model.tes_dispatch_to_pc_con = pe.Constraint(self.model.T,rule=tes_dispatch_to_pc_rule)
        self.model.tes_dispatch_persist_con = pe.Constraint(self.model.T,rule=tes_dispatch_persist_rule)


    def addTESEnergyBalanceConstraints(self):
        """ Method to add TES constraints to the Pyomo Nuclear Model
        
        This method adds constraints pertaining to TES energy balance from charging
        with thermal power and discharging to the power cycle. 
        
        TODO: revisit tes_start_up_rule -> do we need this for nuclear?
        TODO: do we need maintain_tes_rule?
        """
        
        super().addTESEnergyBalanceConstraints()
        
        def tes_balance_rule(model, t):
            """ Balance of energy to and from TES """
            if t == 1:
                return model.s[t] - model.s0 == model.Delta[t] * (model.xntes[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.xtesp[t] + model.Qnsb*model.ynsb[t]))
            return model.s[t] - model.s[t-1] == model.Delta[t] * (model.xntes[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.xtesp[t] + model.Qnsb*model.ynsb[t]))
        def tes_start_up_rule(model, t):
            """ Ensuring sufficient TES charge level to startup NP """
            if t == 1:
                return model.s0 >= model.Delta[t]*model.delta_ns[t]*( (model.Qu + model.Qb)*( -3 + model.ynsu[t] + model.y0 + model.y[t] + model.ycsb0 + model.ycsb[t] ) + model.xtesp[t] + model.Qb*model.ycsb[t] )
            return model.s[t-1] >= model.Delta[t]*model.delta_ns[t]*( (model.Qu + model.Qb)*( -3 + model.ynsu[t] + model.y[t-1] + model.y[t] + model.ycsb[t-1] + model.ycsb[t] ) + model.xtesp[t] + model.Qb*model.ycsb[t] )
        
        # overloaded constraints
        self.model.del_component( self.model.tes_balance_con )
        self.model.del_component( self.model.tes_start_up_con )
        
        self.model.tes_balance_con = pe.Constraint(self.model.T,rule=tes_balance_rule)
        self.model.tes_start_up_con = pe.Constraint(self.model.T,rule=tes_start_up_rule)


    def addNuclearSupplyAndDemandConstraints(self):
        """ Method to add nuclear supply and demand constraints to the Pyomo Nuclear Model
        
        This method adds constraints pertaining to nuclear supply and demand energy
        constraints. Some constraints might be redundant, they are adapted from the CSP
        constraints (thanks LORE team).
        """
        
        super().addNuclearSupplyAndDemandConstraints()
        
        def nuc_production_rule(model,t):
            """ Upper bound on thermal energy produced by NP """
            return model.xnp[t] + model.xntes[t] + model.xnsu[t] + model.Qnsd*model.ynsd[t] <= model.Qin_nuc[t]
        def nuc_generation_rule(model,t):
            """ Thermal energy production by NP only when operating """
            return model.xnp[t] + model.xntes[t] <= model.Qin_nuc[t] * model.yn[t]
        def nuc_min_generation_rule(model,t):
            """ Lower bound on thermal energy produced by NP """
            return model.xnp[t] + model.xntes[t] >= model.Qnl * model.yn[t]
        def nuc_generation_tes_rule(model,t):
            """ Thermal energy production by NP going to TES only when necessary"""
            return model.xntes[t] <= model.Qnhx  * model.yntes[t]
        def nuc_tes_charging_rule(model,t):
            """ Charging of TES via nuclear only when NP operating"""
            return model.yntes[t] <= model.yn[t]
        def nuc_tes_charging_able_rule(model,t):
            """ Charging of TES via nuclear only when NP operating"""
            return model.yntes[t] + model.ytesp[t] <= 1

        # overloaded constraints
        self.model.del_component( self.model.nuc_production_con )
        self.model.del_component( self.model.nuc_generation_con )
        self.model.del_component( self.model.nuc_min_generation_con )
        
        self.model.nuc_production_con = pe.Constraint(self.model.T,rule=nuc_production_rule)
        self.model.nuc_generation_con = pe.Constraint(self.model.T,rule=nuc_generation_rule)
        self.model.nuc_min_generation_con = pe.Constraint(self.model.T,rule=nuc_min_generation_rule)
        
        # new constraints
        self.model.nuc_generation_tes_con = pe.Constraint(self.model.T,rule=nuc_generation_tes_rule)
        self.model.nuc_tes_charging_con = pe.Constraint(self.model.T,rule=nuc_tes_charging_rule)
        self.model.nuc_tes_charging_able_rule_con = pe.Constraint(self.model.T,rule=nuc_tes_charging_able_rule)


    def generate_constraints(self, skip_parent=False):
        """ Method to add ALL constraints to the Pyomo Nuclear Model
        
        This method calls the previously defined constraint methods to instantiate
        them and add to the existing model. This method first calls the GeneralDispatch
        version to set PowerCycle constraints, then calls nuclear constraint methods
        to add them to the model. 
        """
        
        super().generate_constraints()

        
# =============================================================================
# Dispatch Wrapper
# =============================================================================

class IndirectNuclearDispatchParamWrap(NuclearDispatchParamWrap):
    """
    The NuclearDispatchParamWrap class is meant to be the staging area for the 
    creation of Parameters ONLY for the NuclearDispatch class. It communicates 
    with the NE2 modules, receiving SSC and PySAM input dictionaries to calculate 
    both static parameters used for every simulation segment AND initial conditions 
    that can be updated.
    """
    
    def __init__(self, **kwargs):
        """ Initializes the NuclearDispatchParamWrap module
        
        Inputs:
            unitRegistry (pint.registry)   : unique unit Pint unit registry
            SSC_dict (dict)                : dictionary of SSC inputs needed to run modules
            PySAM_dict (dict)              : dictionary of PySAM inputs + file names
            pyomo_horizon (int Quant)      : length of Pyomo simulation segment (hours)
            dispatch_time_step (int Quant) : length of each Pyomo time step (hours)
        """
        kwargs['direct'] = kwargs['direct'] if 'direct' in kwargs else False
        super().__init__(**kwargs)


    def set_design(self):
        """ Method to calculate and save design point values of Plant operation
        
        This method extracts values and calculates for design point parameters 
        of our Plant (e.g., nuclear thermal power output, power cycle efficiency,
        inlet and outlet temperatures, etc.). 
        """
        
        super(NuclearDispatchParamWrap, self ).set_design()
        
        u = self.u
        
        # nuclear parameters
        self.q_nuc_design = self.SSC_dict['q_dot_nuclear_des'] * u.MW      # nuclear design thermal power
        
        # specific heat values at design point
        T_htf  = 0.5*(self.T_htf_hot + self.T_htf_cold)
        cp_des = SSCHelperMethods.get_cp_htf(self.u, T_htf, self.SSC_dict['nuc_htf'], self.interpolants['cp_interp'] )
        cp_des = cp_des.to('J/g/kelvin')       
        
        # mass flow rate
        dm_des = self.q_pb_design / (cp_des * (self.T_htf_hot - self.T_htf_cold) )  
        self.dm_pb_design = dm_des.to('kg/s')                               # power block design mass flow rate
        
        # TES design point
        e_tes_design = self.q_pb_design * self.SSC_dict['tshours']*u.hr  
        m_tes_des = e_tes_design / cp_des / (self.T_htf_hot - self.T_htf_cold)     
        self.e_tes_design = e_tes_design.to('kWh') # TES storage capacity (kWht)
        self.m_tes_design = m_tes_des.to('kg')     # TES active storage mass (kg)


    def set_indirect_config_parameters(self, param_dict):
        """ Method to set fixed costs of the Plant
        
        This method calculates some fixed costs for the Plant operations, startup,
        standby, etc. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        u = self.u

        Qnc     = self.q_nuc_design
        propn_LFR_flow = (self.q_nuc_design / self.q_pb_design).m
        slope,intercept = SSCHelperMethods.linearize_indirectTES_eff(propn_LFR_flow)
        
        eta_lin  = slope * self.eta_design
        intcpts  = intercept * self.p_pb_design
        
        Wnc     = Qnc * self.eta_design  
        Qnhx    = self.q_nuc_design * 0.9  # TODO: simulating hit in efficiency for nuclear thermal power when charging TES

        ### Cost Parameters ###
        param_dict['Wnc']    = Wnc.to('kW')       #W^{nc}: Cycle capacity for nuclear power only [kWe]
        param_dict['Qnhx']   = Qnhx.to('kW')      #Q^{nhx}: Upper limit on nuclear power to charge TES due to HX efficiency [kWt]
        param_dict['eta_LD'] = eta_lin[1]         #\eta^{LD}: Linearized cycle efficiency during low demand operation [-]
        param_dict['eta_HD'] = eta_lin[0]         #\eta^{HD}: Linearized cycle efficiency during high demand operation [-]
        param_dict['intcpts_LD'] = intcpts[1].to('kW')         #\eta^{LD}: Linearized cycle efficiency during low demand operation [-]
        param_dict['intcpts_HD'] = intcpts[0].to('kW')         #\eta^{HD}: Linearized cycle efficiency during high demand operation [-]
        
        return param_dict

# =============================================================================
# Dispatch Outputs
# =============================================================================
  
class IndirectNuclearDispatchOutputs(object):
    """
    The IndirectNuclearDispatchOutputs class is meant to handle outputs from a given,
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
        Qnhx = pe.value(dm.model.Qnhx)/1000.
        # nuclear power going to PC turbine
        xn = np.array([pe.value(dm.model.xnp[t])   for t in t_pyomo])/1000. # from kWt -> MWt
        # total power into PC (nuclear + TES dispatch)
        x = xn + np.array([pe.value(dm.model.xtesp[t])   for t in t_pyomo])/1000. # from kWt -> MWt
        # nuclear power going to TES charging
        xntes = np.array([pe.value(dm.model.xntes[t])   for t in t_pyomo])/1000.
        
        #----Thermal Capacity for Cycle Startup and Operation----
        Qc = np.array([pe.value(dm.model.Qc[t]) for t in t_pyomo])/1000. # from kWt -> MWt
        Qn = np.array([pe.value(dm.model.Qin_nuc[t]) for t in t_pyomo])/1000. # from kWt -> MWt
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
        f_nuc_to_tes_target  = [xntes[t] / Qnhx if xn[t] < 0.001 else xntes[t] / Qn[t] for t in t_horizon]
        
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
            disp_targs['f_nuc_to_tes_target']  = f_nuc_to_tes_target
            
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
            disp_targs['f_nuc_to_tes_target']  = np.hstack( [f_nuc_to_tes_target  , empty_array] ).tolist()
            
        return disp_targs