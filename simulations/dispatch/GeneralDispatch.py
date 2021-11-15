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
import os, copy
import pyomo.environ as pe
import numpy as np
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
from pyomo.environ import units as u_pyomo
if not hasattr(u_pyomo,'USD'):
    u_pyomo.load_definitions_from_strings(['USD = [currency]'])
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
from abc import ABC, abstractmethod


class GeneralDispatch(ABC):
    """
    The GeneralDispatch class is meant to set up and run Dispatch
    optimization as a mixed integer linear program problem using Pyomo.
    It creates a ConcreteModel in Pyomo with Parameters, Variables,
    Objectives, and Constraints. The ConcreteModel can then be solved. 
    
    This is the base/parent class and should not be run by itself. 
    Essentially, this sets Pyomo with PowerCycle-specific parameters,
    variables, objectives, and constraints. Derived classes can then
    be instantiated to add more specific constraints and overload the
    objective function. 
    
    TODO: can we convert this to an abstract class in Python?
    """
    
    @abstractmethod
    def __init__(self, params, unitRegistry):
        """ Initializes the GeneralDispatch module
        
        The instantiation of this class receives a parameter dictionary from
        the NE2 module (created using a DispatchWrapper class in this file but
        separate from the GeneralDispatch or other Dispatch classes). It first
        creates an empty Concrete Model from Pyomo, then generates Parameters
        from the parameter dictionary, Variables, Objectives and Constraints.
        
        Inputs:
            params (dict)                : dictionary of Pyomo dispatch parameters
            unitRegistry (pint.registry) : unique unit Pint unit registry
        """
        
        self.u_pyomo = u_pyomo
        
        self.model = pe.ConcreteModel()
        self.generate_params(params)
        self.generate_variables()
        self.add_objective()
        self.generate_constraints()
        
        # assert_units_consistent(self.model)


    def generate_params(self, params):
        """ Method to generate parameters within Pyomo General Model
        
        This method reads in a dictionary of pyomo parameters and uses these
        inputs to initialize parameters for the Pyomo Concrete Model. This method
        sets up parameters particularly for the Power Cycle. It also defines
        some lambda functions that helps convert Pint units to Pyomo units.
        
        Note: initial conditions are defined for the time period immediately 
        preceding the start of this new Pyomo time segment. 
        
        Inputs:
            params (dict)  : dictionary of Pyomo dispatch parameters
        """
        
        time_range = range(1,params["T"]+1)
        # this is a lambda function to grab Pint Quantity magntitude if it has one (gm = get magnitude)
        self.gm = lambda param_str: params[param_str].m if hasattr(params[param_str],'m') else params[param_str]
        gm = self.gm
        # creates dictionary if parameter happens to be a time-indexed parameter (gd = get dictionary)
        self.gd = lambda param_str: {ind:val for (ind,val) in zip(time_range, gm(param_str))} if type(gm(param_str)) is np.ndarray else gm(param_str)
        gd = self.gd
        # shortening definition name for converting Pint units to pyomo units (gu = get unit)
        self.gu = lambda param_str: SSCHelperMethods.convert_to_pyomo_unit(params, param_str)
        gu = self.gu
        
        ### Time Indexed Parameters ###
        self.model.T = pe.Set(initialize = time_range)            #T: time periods
        self.model.num_periods = pe.Param(initialize=params["T"]) #N_T: number of time periods
        self.model.Delta = pe.Param(self.model.T, mutable=False, initialize=gd("Delta"), units=gu("Delta"))          #\Delta_{t}: duration of period t
        self.model.Delta_e = pe.Param(self.model.T, mutable=False, initialize=gd("Delta_e"), units=gu("Delta_e"))    #\Delta_{e,t}: cumulative time elapsed at end of period t
        self.model.D = pe.Param(self.model.T, mutable=True, initialize=gd("D"), units=gu("D"))                          #D_{t}: Time-weighted discount factor in period $t$ [-]
        self.model.etaamb = pe.Param(self.model.T, mutable=True, initialize=gd("etaamb"), units=gu("etaamb"))           #\eta^{amb}_{t}: Cycle efficiency ambient temperature adjustment factor in period $t$ [-]
        self.model.etac = pe.Param(self.model.T, mutable=True, initialize=gd("etac"), units=gu("etac"))                 #\eta^{c}_{t}: Normalized condenser parasitic loss in period $t$ [-] 
        self.model.P = pe.Param(self.model.T, mutable=True, initialize=gd("P"), units=gu("P"))                          #P_{t}: Electricity sales price in period $t$ [\$/kWh]
        self.model.Qc = pe.Param(self.model.T, mutable=True, initialize=gd("Qc"), units=gu("Qc"))                       #Q^{c}_{t}: Allowable power per period for cycle start-up in period $t$ [kWt]
        self.model.Wdotnet = pe.Param(self.model.T, mutable=True, initialize=gd("Wdotnet"), units=gu("Wdotnet"))        #\dot{W}^{net}_{t}: Net grid transmission upper limit in period $t$ [kWe]
        self.model.W_u_plus = pe.Param(self.model.T, mutable=True, initialize=gd("W_u_plus"), units=gu("W_u_plus"))     #W^{u+}_{t}: Maximum power production when starting generation in period $t$  [kWe]
        self.model.W_u_minus = pe.Param(self.model.T, mutable=True, initialize=gd("W_u_minus"), units=gu("W_u_minus"))  #W^{u-}_{t}: Maximum power production in period $t$ when stopping generation in period $t+1$  [kWe]

        ### Power Cycle Parameters ###
        self.model.Ec = pe.Param(mutable=True, initialize=gd("Ec"), units=gu("Ec"))           #E^c: Required energy expended to start cycle [kWt$\cdot$h]
        self.model.eta_des = pe.Param(mutable=True, initialize=gd("eta_des"), units=gu("eta_des"))   #\eta^{des}: Cycle nominal efficiency [-] 
        self.model.etap = pe.Param(mutable=True, initialize=gd("etap"), units=gu("etap"))     #\eta^p: Slope of linear approximation of power cycle performance curve [kWe/kWt]
        self.model.Lc = pe.Param(mutable=True, initialize=gd("Lc"), units=gu("Lc"))           #L^c: Cycle heat transfer fluid pumping power per unit energy expended [kWe/kWt]
        self.model.Qb = pe.Param(mutable=True, initialize=gd("Qb"), units=gu("Qb"))           #Q^b: Cycle standby thermal power consumption per period [kWt]
        self.model.Ql = pe.Param(mutable=True, initialize=gd("Ql"), units=gu("Ql"))           #Q^l: Minimum operational thermal power input to cycle [kWt]
        self.model.Qu = pe.Param(mutable=True, initialize=gd("Qu"), units=gu("Qu"))           #Q^u: Cycle thermal power capacity [kWt]
        self.model.Wb = pe.Param(mutable=True, initialize=gd("Wb"), units=gu("Wb"))           #W^b: Power cycle standby operation parasitic load [kWe]
        self.model.Wdotl = pe.Param(mutable=True, initialize=gd("Wdotl"), units=gu("Wdotl"))          #\dot{W}^l: Minimum cycle electric power output [kWe]
        self.model.Wdotu = pe.Param(mutable=True, initialize=gd("Wdotu"), units=gu("Wdotu"))          #\dot{W}^u: Cycle electric power rated capacity [kWe]
        self.model.W_delta_plus = pe.Param(mutable=True, initialize=gd("W_delta_plus"), units=gu("W_delta_plus"))    #W^{\Delta+}: Power cycle ramp-up designed limit [kWe/h]
        self.model.W_delta_minus = pe.Param(mutable=True, initialize=gd("W_delta_minus"), units=gu("W_delta_minus")) #W^{\Delta-}: Power cycle ramp-down designed limit [kWe/h]
        self.model.W_v_plus = pe.Param(mutable=True, initialize=gd("W_v_plus"), units=gu("W_v_plus"))       #W^{v+}: Power cycle ramp-up violation limit [kWe/h]
        self.model.W_v_minus = pe.Param(mutable=True, initialize=gd("W_v_minus"), units=gu("W_v_minus"))    #W^{v-}: Power cycle ramp-down violation limit [kWe/h]
        self.model.Yu = pe.Param(mutable=True, initialize=gd("Yu"), units=gu("Yu"))           #Y^u: Minimum required power cycle uptime [h]
        self.model.Yd = pe.Param(mutable=True, initialize=gd("Yd"), units=gu("Yd"))           #Y^d: Minimum required power cycle downtime [h]

        ### Cost Parameters ###
        self.model.alpha = pe.Param(mutable=True, initialize=gd("alpha"), units=gu("alpha"))           #\alpha: Conversion factor between unitless and monetary values [\$]
        self.model.Cpc = pe.Param(mutable=True, initialize=gd("Cpc"), units=gu("Cpc"))                 #C^{pc}: Operating cost of power cycle [\$/kWe$\cdot$h]
        self.model.Ccsu = pe.Param(mutable=True, initialize=gd("Ccsu"), units=gu("Ccsu"))              #C^{csu}: Penalty for power cycle cold start-up [\$/start]
        self.model.Cchsp = pe.Param(mutable=True, initialize=gd("Cchsp"), units=gu("Cchsp"))           #C^{chsp}: Penalty for power cycle hot start-up [\$/start]
        self.model.C_delta_w = pe.Param(mutable=True, initialize=gd("C_delta_w"), units=gu("C_delta_w"))    #C^{\delta_w}: Penalty for change in power cycle  production [\$/$\Delta\text{kWe}$]
        self.model.C_v_w = pe.Param(mutable=True, initialize=gd("C_v_w"), units=gu("C_v_w"))           #C^{v_w}: Penalty for change in power cycle  production tcb{beyond designed limits} [\$/$\Delta\text{kWe}$]
        self.model.Ccsb = pe.Param(mutable=True, initialize=gd("Ccsb"), units=gu("Ccsb"))              #C^{csb}: Operating cost of power cycle standby operation [\$/kWt$\cdot$h]
        
        ### Initial Condition Parameters ###
        self.model.ucsu0 = pe.Param(mutable=True, initialize=gd("ucsu0"), units=gu("ucsu0")) #u^{csu}_0: Initial cycle start-up energy inventory  [kWt$\cdot$h]
        self.model.wdot0 = pe.Param(mutable=True, initialize=gd("wdot0"), units=gu("wdot0")) #\dot{w}_0: Initial power cycle electricity generation [kW]e
        self.model.y0 = pe.Param(mutable=True, initialize=gd("y0"), units=gu("y0"))          #y_0: 1 if cycle is generating electric power initially, 0 otherwise
        self.model.ycsb0 = pe.Param(mutable=True, initialize=gd("ycsb0"), units=gu("ycsb0")) #y^{csb}_0: 1 if cycle is in standby mode initially, 0 otherwise
        self.model.ycsu0 = pe.Param(mutable=True, initialize=gd("ycsu0"), units=gu("ycsu0")) #y^{csu}_0: 1 if cycle is in starting up initially, 0 otherwise    [az] this is new.
        self.model.Yu0 = pe.Param(mutable=True, initialize=gd("Yu0"), units=gu("Yu0"))       #Y^u_0: duration that cycle has been generating electric power [h]
        self.model.Yd0 = pe.Param(mutable=True, initialize=gd("Yd0"), units=gu("Yd0"))       #Y^d_0: duration that cycle has not been generating power (i.e., shut down or in standby mode) [h]
        
        #------- Persistence Parameters ---------
        # TODO: removing references to wdot_s_prev, they only exist as Constraints and should later be added to Objective somehow
        # self.model.wdot_s_prev  = pe.Param(self.model.T, mutable=True, initialize=gd("wdot_s_prev"), units=gu("wdot_s_prev")) #\dot{w}^{s,prev}: previous $\dot{w}^s$, or energy sold to grid [kWe]
        # self.model.wdot_s_pen  = pe.Param(self.model.T, mutable=True, initialize=gd("wdot_s_pen"), units=gu("delta_rs"))    #\dot{w}_{s,pen}: previous $\dot{w}$ 


    def generate_variables(self):
        """ Method to generate parameters within Pyomo General Model
        
        This method instantiates variables for the Pyomo Concrete Model, with
        domains. Does not need initial guesses here, they are defined in the 
        parameters. Here we define continuous and binary variables for the 
        Power Cycle. 
        """
        u = self.u_pyomo
        ### Decision Variables ###
        #------- Variables ---------
        self.model.ucsu = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kWh)                             #u^{csu}: Cycle start-up energy inventory at period $t$ [kWt$\cdot$h]
        self.model.wdot = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)                              #\dot{w}: Power cycle electricity generation at period $t$ [kWe]
        self.model.wdot_delta_plus = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW/u.hr)	                #\dot{w}^{\Delta+}: Power cycle ramp-up in period $t$ [kWe/hr]
        self.model.wdot_delta_minus = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW/u.hr) 	                #\dot{w}^{\Delta-}: Power cycle ramp-down in period $t$ [kWe/hr]
        self.model.wdot_v_plus = pe.Var(self.model.T, domain=pe.NonNegativeReals, bounds = (0,self.model.W_v_plus), units=u.kW/u.hr)      #\dot{w}^{v+}: Power cycle ramp-up beyond designed limit in period $t$ [kWe/hr]
        self.model.wdot_v_minus = pe.Var(self.model.T, domain=pe.NonNegativeReals, bounds = (0,self.model.W_v_minus), units=u.kW/u.hr) 	 #\dot{w}^{v-}: Power cycle ramp-down beyond designed limit in period $t$ [kWe/hr]
        self.model.wdot_s = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)	                       #\dot{w}^s: Energy sold to grid in time t [kWe]
        self.model.wdot_p = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)	                       #\dot{w}^p: Energy purchased from the grid in time t [kWe]
        self.model.x = pe.Var(self.model.T, domain=pe.NonNegativeReals, units=u.kW)                                #x: Cycle thermal power utilization at period $t$ [kWt]
        
        #------- Binary Variables ---------
        self.model.y = pe.Var(self.model.T, domain=pe.Binary)         #y: 1 if cycle is generating electric power at period $t$; 0 otherwise
        self.model.ychsp = pe.Var(self.model.T, domain=pe.Binary)     #y^{chsp}: 1 if cycle hot start-up penalty is incurred at period $t$ (from standby); 0 otherwise
        self.model.ycsb = pe.Var(self.model.T, domain=pe.Binary)      #y^{csb}: 1 if cycle is in standby mode at period $t$; 0 otherwise
        self.model.ycsd = pe.Var(self.model.T, domain=pe.Binary)	  #y^{csd}: 1 if cycle is shutting down at period $t$; 0 otherwise
        self.model.ycsu = pe.Var(self.model.T, domain=pe.Binary)      #y^{csu}: 1 if cycle is starting up at period $t$; 0 otherwise
        self.model.ycsup = pe.Var(self.model.T, domain=pe.Binary)     #y^{csup}: 1 if cycle cold start-up penalty is incurred at period $t$ (from off); 0 otherwise
        self.model.ycgb = pe.Var(self.model.T, domain=pe.NonNegativeReals, bounds=(0,1))      #y^{cgb}: 1 if cycle begins electric power generation at period $t$; 0 otherwise
        self.model.ycge = pe.Var(self.model.T, domain=pe.NonNegativeReals, bounds=(0,1))      #y^{cge}: 1 if cycle stops electric power generation at period $t$; 0 otherwise
        
        #------- Persistence Variables ---------
        # TODO: revisit wdot_s_prev references later
        # self.model.wdot_s_prev_delta_plus = pe.Var(self.model.T, domain=pe.NonNegativeReals)  #\dot{w}^{\Delta+}_{s,prev}: upper bound on energy sold [kWe]
        # self.model.wdot_s_prev_delta_minus = pe.Var(self.model.T, domain=pe.NonNegativeReals) #\dot{w}^{\Delta-}_{s,prev}: lower bound on energy sold [kWe]
        # self.model.ycoff = pe.Var(self.model.T, domain=pe.Binary)     #y^{c,off}: 1 if power cycle is off at period $t$; 0 otherwise
        
    @abstractmethod
    def add_objective(self):
        """ Method to add an objective function to the Pyomo General Model
        
        This method adds an objective function to the Pyomo General Dispatch
        model. Typically, the LORE team defined a nested function and passes it
        into the Pyomo model. This version of the method should not be called,
        is designated as an abstract method.
        """
        
        def objectiveRule(model):
            """ Maximize profits """
            return (
                    sum( model.D[t] * model.Delta[t] * model.P[t]
                    for t in model.T) 
                    )
        
        self.model.OBJ = pe.Objective(rule=objectiveRule, sense = pe.maximize)


    def addPersistenceConstraints(self):
        """ Method to add persistence constraints to the Pyomo General Model
        
        This method adds constraints pertaining to persistence within the Pyomo
        General Dispatch class. Namely, it constraints the amount of energy
        sold into the grid in timestep t+1 relative to timestep t. It makes sure
        that amount doesn't fluctuate by too much, whether increansig or decreasing.
        """
        
        # TODO: wdot_s_prev should be a single float, there should be an {if t==1} statement for wdot_s_prev
        #              and for all other t, we compare wdot_s[t] - wdot_s[t-1]
        def wdot_s_persist_pos_rule(model,t):
            """ Energy sold shouldn't increase too much """
            return model.wdot_s_prev_delta_plus[t] >= model.wdot_s[t] - model.wdot_s_prev[t]
        def wdot_s_persist_neg_rule(model,t):
            """ Energy sold shouldn't decrease too much """
            return model.wdot_s_prev_delta_minus[t] >= model.wdot_s_prev[t] - model.wdot_s[t]
        
        self.model.persist_pos_con = pe.Constraint(self.model.T,rule=wdot_s_persist_pos_rule)
        self.model.persist_neg_con = pe.Constraint(self.model.T,rule=wdot_s_persist_neg_rule)
        
        
    def addCycleStartupConstraints(self):
        """ Method to add cycle startup constraints to the Pyomo General Model
        
        This method adds constraints pertaining to cycle startup within the Pyomo
        General Dispatch class. Several nested functions are defined. 
        """
        def pc_inventory_rule(model, t):
            """ PC energy inventory for startup is balanced """
            if t == 1:
                return model.ucsu[t] <= model.ucsu0 + model.Delta[t] * model.Qc[t] * model.ycsu[t]
            return model.ucsu[t] <= model.ucsu[t-1] + model.Delta[t] * model.Qc[t] * model.ycsu[t]
        def pc_inv_nonzero_rule(model, t):
            """ PC energy inventory is positive if cycle starting up """
            return model.ucsu[t] <= model.Ec * model.ycsu[t]
        def pc_startup_rule(model, t):
            """ PC resumes normal operation after startup, standby, or previously operating """
            current_Delta = model.Delta[t]
            # structure of inequality is lb <= model.param <= ub with strict=False by default
            if self.eval_ineq(1,current_Delta) and t == 1:
                return model.y[t] <= model.ucsu[t]/model.Ec + model.y0 + model.ycsb0
            elif self.eval_ineq(1,current_Delta) and t > 1:
                return model.y[t] <= model.ucsu[t]/model.Ec + model.y[t-1] + model.ycsb[t-1]
            elif self.eval_ineq(current_Delta,1,strict=True) and t == 1:
                return self.eval_ineq(model.y[t], model.ucsu0/model.Ec + model.y0 + model.ycsb0)
            # only case remaining: Delta[t]<1, t>1
            return self.eval_ineq(model.y[t], model.ucsu[t-1]/model.Ec + model.y[t-1] + model.ycsb[t-1])
        def pc_production_rule(model, t):
            """ Heat input used for PC startup doesn't exceed max """
            return model.x[t] + model.Qc[t]*model.ycsu[t] <= model.Qu
        def pc_generation_rule(model, t):
            """ Upper bound on heat input to PC """
            return model.x[t] <= model.Qu * model.y[t]
        def pc_min_gen_rule(model, t):
            """ Lower bound on heat input to PC """
            return model.x[t] >= model.Ql * model.y[t]
        
        self.model.pc_inventory_con = pe.Constraint(self.model.T,rule=pc_inventory_rule)
        self.model.pc_inv_nonzero_con = pe.Constraint(self.model.T,rule=pc_inv_nonzero_rule)
        self.model.pc_startup_con = pe.Constraint(self.model.T,rule=pc_startup_rule)
        self.model.pc_production_con = pe.Constraint(self.model.T,rule=pc_production_rule)
        self.model.pc_generation_con = pe.Constraint(self.model.T,rule=pc_generation_rule)
        self.model.pc_min_gen_con = pe.Constraint(self.model.T,rule=pc_min_gen_rule)
        
        
    def addPiecewiseLinearEfficiencyConstraints(self):
        """ Method to add efficiency constraints to the Pyomo General Model
        
        This method adds constraints pertaining to efficiency constraints defined
        as a piecewise linear approximation. Also referred to as Cycle supply and 
        demand constraints.
        """
        def power_rule(model, t):
            """ Model of electric power vs. heat input as linear function """
            return model.wdot[t] == (model.etaamb[t]/model.eta_des)*(model.etap*model.x[t] + model.y[t]*(model.Wdotu - model.etap*model.Qu))
        def power_ub_rule(model, t):
            """ Upper bound on PC electric power output """
            return model.wdot[t] <= model.Wdotu*(model.etaamb[t]/model.eta_des)*model.y[t]
        def power_lb_rule(model, t):
            """ Lower bound on PC electric power output """
            return model.wdot[t] >= model.Wdotl*(model.etaamb[t]/model.eta_des)*model.y[t]
        def change_in_w_pos_rule(model, t):
            """ Limit on increase in electric power production """
            if t == 1:
                return model.wdot_delta_plus[t] >= model.wdot[t] - model.wdot0
            return model.wdot_delta_plus[t] >= model.wdot[t] - model.wdot[t-1]
        def change_in_w_neg_rule(model, t):
            """ Limit on decrease in electric power production """
            if t == 1:
                return model.wdot_delta_minus[t] >= model.wdot0 - model.wdot[t]
            return model.wdot_delta_minus[t] >= model.wdot[t-1] - model.wdot[t]
        def cycle_ramp_rate_pos_rule(model, t):
            """ Limit on excessive ramp up of PC power production """
            return (
                    model.wdot_delta_plus[t] - model.wdot_v_plus[t] <= model.W_delta_plus*model.Delta[t] 
                    + ((model.etaamb[t]/model.eta_des)*model.W_u_plus[t] - model.W_delta_plus*model.Delta[t])*model.ycgb[t]
            )
        def cycle_ramp_rate_neg_rule(model, t):
            """ Limit on excessive ramp down of PC power production"""
            return (
                    model.wdot_delta_minus[t] - model.wdot_v_minus[t] <= model.W_delta_minus*model.Delta[t] 
                    + ((model.etaamb[t]/model.eta_des)*model.W_u_minus[t] - model.W_delta_minus*model.Delta[t])*model.ycge[t]
            )
        def grid_max_rule(model, t):
            """ Limits grid transmission for net power production """
            return model.wdot_s[t] <= model.Wdotnet[t]

        
        self.model.power_con = pe.Constraint(self.model.T,rule=power_rule)
        self.model.power_ub_con = pe.Constraint(self.model.T,rule=power_ub_rule)
        self.model.power_lb_con = pe.Constraint(self.model.T,rule=power_lb_rule)
        self.model.change_in_w_pos_con = pe.Constraint(self.model.T,rule=change_in_w_pos_rule)
        self.model.change_in_w_neg_con = pe.Constraint(self.model.T,rule=change_in_w_neg_rule)
        self.model.cycle_ramp_rate_pos_con = pe.Constraint(self.model.T,rule=cycle_ramp_rate_pos_rule)
        self.model.cycle_ramp_rate_neg_con = pe.Constraint(self.model.T,rule=cycle_ramp_rate_neg_rule)
        self.model.grid_max_con = pe.Constraint(self.model.T,rule=grid_max_rule)
        
        
    def addMinUpAndDowntimeConstraints(self):
        """ Method to add up and downtime constraints to the Pyomo General Model
        
        This method adds constraints pertaining to enforcing minimum uptime
        and downtime of the power cycle whenever it starts or ends operation,
        respectively. 
        """
        def min_cycle_uptime_rule(model,t):
            """ Enforce minimum cycle uptime """
            if pe.value(model.Delta_e[t] > (model.Yu - model.Yu0) * model.y0):
                return sum(model.ycgb[tp] for tp in model.T if pe.value(model.Delta_e[t]-model.Delta_e[tp] < model.Yu) and pe.value(model.Delta_e[t] - model.Delta_e[tp] >= 0)) <= model.y[t]
            return pe.Constraint.Feasible
        def min_cycle_downtime_rule(model,t):
            """ Enforce minimum cycle downtime """
            if pe.value(model.Delta_e[t] > ((model.Yd - model.Yd0)*(1-model.y0))):
                return sum( model.ycge[tp] for tp in model.T if pe.value(model.Delta_e[t]-model.Delta_e[tp] < model.Yd) and pe.value(model.Delta_e[t] - model.Delta_e[tp] >= 0))  <= (1 - model.y[t])
            return pe.Constraint.Feasible
        def cycle_start_end_gen_rule(model,t):
            """ Tracks startup and shutdown of plants (Garver 1962) """
            if t == 1:
                return model.ycgb[t] - model.ycge[t] == model.y[t] - model.y0
            return model.ycgb[t] - model.ycge[t] == model.y[t] - model.y[t-1]
        def cycle_min_updown_init_rule(model,t):
            """ Initial conditions for minimum uptime/downtime constraints """
            if self.eval_ineq(model.Delta_e[t],
                              max(pe.value(model.y0*(model.Yu-model.Yu0)), pe.value((1-model.y0)*(model.Yd-model.Yd0)))):
                return model.y[t] == model.y0
            return pe.Constraint.Feasible
        
        self.model.min_cycle_uptime_con = pe.Constraint(self.model.T,rule=min_cycle_uptime_rule)
        self.model.min_cycle_downtime_con = pe.Constraint(self.model.T,rule=min_cycle_downtime_rule)
        self.model.cycle_start_end_gen_con = pe.Constraint(self.model.T,rule=cycle_start_end_gen_rule)
        self.model.cycle_min_updown_init_con = pe.Constraint(self.model.T,rule=cycle_min_updown_init_rule)
        
        
    def addCycleLogicConstraints(self):
        """ Method to add cycle logic constraints to the Pyomo General Model
        
        This method adds constraints pertaining to tracking binary variable
        logic for the power cycle. Essentially, to make sure the correct modes
        are activated when they are allowed. 
        """
        def pc_su_persist_rule(model, t):
            """ PC startup can't follow an operating timestep """
            if t == 1:
                return model.ycsu[t] + model.y0 <= 1
            return model.ycsu[t] + model.y[t-1] <= 1
        def pc_su_subhourly_rule(model, t):
            """ PC startup and operation can't coincide """
            if self.eval_ineq(model.Delta[t],1,strict=True):
                return model.y[t] + model.ycsu[t] <= 1
            return pe.Constraint.Feasible  #no analogous constraint for hourly or longer time steps
        def pc_sb_start_rule(model, t):
            """ PC standby must follow operating or another standby timestep """
            if t == 1:
                return model.ycsb[t] <= model.y0 + model.ycsb0
            return model.ycsb[t] <= model.y[t-1] + model.ycsb[t-1]
        def pc_sb_part1_rule(model, t):
            """ PC startup and standby can't coincide """
            return model.ycsu[t] + model.ycsb[t] <= 1
        def pc_sb_part2_rule(model, t):
            """ PC operation and standby can't coincide """
            return model.y[t] + model.ycsb[t] <= 1
        def cycle_sb_pen_rule(model, t):
            """ PC hot startup penalty: if we went from standby to operating """
            if t == 1:
                 return model.ychsp[t] >= model.y[t] - (1 - model.ycsb0)
            return model.ychsp[t] >= model.y[t] - (1 - model.ycsb[t-1])
        def cycle_shutdown_rule(model, t):
            """ PC shutdown only after operating or standby """
            if t == 1:
                return model.ycsd[t] >= model.y0 - model.y[t] + model.ycsb0 - model.ycsb[t]
            return model.ycsd[t] >= model.y[t-1] - model.y[t] + model.ycsb[t-1] - model.ycsb[t]
        def cycle_start_pen_rule(model, t):
            """ PC cold startup penalty """
            if t == 1: 
                return model.ycsup[t] >= model.ycsu[t] - model.ycsu0 
            return model.ycsup[t] >= model.ycsu[t] - model.ycsu[t-1]
         
        self.model.pc_su_persist_con = pe.Constraint(self.model.T,rule=pc_su_persist_rule)
        self.model.pc_su_subhourly_con = pe.Constraint(self.model.T,rule=pc_su_subhourly_rule)
        self.model.pc_sb_start_con = pe.Constraint(self.model.T,rule=pc_sb_start_rule)
        self.model.pc_sb_part1_con = pe.Constraint(self.model.T,rule=pc_sb_part1_rule)
        self.model.pc_sb_part2_con = pe.Constraint(self.model.T,rule=pc_sb_part2_rule)
        self.model.cycle_sb_pen_con = pe.Constraint(self.model.T,rule=cycle_sb_pen_rule)
        self.model.cycle_shutdown_con = pe.Constraint(self.model.T,rule=cycle_shutdown_rule)
        self.model.cycle_start_pen_con = pe.Constraint(self.model.T,rule=cycle_start_pen_rule)


    def generate_constraints(self):
        """ Method to add ALL constraints to the Pyomo General Model
        
        This method calls the previously defined constraint methods to instantiate
        them and add to the existing model. This method can be called from a 
        derived class as well as overloaded.
        """
        
        # self.addPersistenceConstraints() # TODO: revisit these constraints later when we know how to add a term to Objective
        self.addPiecewiseLinearEfficiencyConstraints()
        self.addCycleStartupConstraints()
        self.addMinUpAndDowntimeConstraints()
        self.addCycleLogicConstraints()
    
    
    def eval_ineq(self, lb, expr, ub=None, strict=False, val=True):
        """ Method to evaluate an inequality using pyomo syntax
        
        This is a helper method to evaluate an inequality using the specific
        pyomo syntax needed to evaluate when solving the model. 
        
        Inputs:
            lb (**)       : lower bound of inequality
            expr (**)     : middle value of inequality
            ub (**)       : upper bound of inequality (if exists)
            strict (bool) : strictness of inequality (True means < or >)
            val (bool)    : flag to return value instead of expression
        Outputs:
            ineq (logical expr or float) : evaluated expression as a Pyomo object or float
            
        (**) - can either be a float input or a ParamData object
        """
        ineq = pe.inequality(lb, expr, ub, strict=strict)
        
        # return value if requested, otherwise return Inequality Expression
        if val:
            return pe.value(ineq)
        else:
            return ineq


    def solve_model(self, mipgap=0.7, tee=False, run_simple=False):
        """ Method to solve the Pyomo model
        
        This method solves the Pyomo Concrete model that has been instantiated
        and constructed in the __init__ of this class. The MILP (mixed integer
        linear programming) problem is solved using a CBC (coin-or branch and cut)
        solver that should be installed prior to solving. It takes in an input
        for the MIP or Ratio gap. It is a termination condition for the optimization:
        if the difference between the best known solution and worst known solution
        (lower and upper bounds, respectively) is lower than this gap, we have
        found our optimal solution. Here is a good primer:
            https://www.gurobi.com/resource/mip-basics/
        
        Inputs:
            mipgap (float): 
                minimum gap to find optimal solution
        Outputs:
            results (SolverResults) : 
                dictionary of solver output after optimization, contains 
                information on convergence, etc.
        """
        
        # define solver (Coin-or branch and cut)
        opt = pe.SolverFactory('cbc')
        
        # setting optimality condition
        opt.options["ratioGap"] = mipgap
        
        if run_simple:
            opt.options["primalPivot"] = "dantzig"
            opt.options["dualPivot"]   = "dantzig"
        
        # solving model
        results = opt.solve(self.model, tee=tee, keepfiles=False)
        #TODO: assert successful optimization?
        
        return results
    
    
# =============================================================================
# Dispatch Wrapper
# =============================================================================

class GeneralDispatchParamWrap(object):
    """
    The GeneralDispatchParamWrap class is meant to be the staging area for the 
    creation of Parameters ONLY. It communicates with the NE2 modules, receiving
    SSC and PySAM input dictionaries to calculate both static parameters used 
    for every simulation segment AND initial conditions that can be updated.
    """
    
    def __init__(self, unit_registry, SSC_dict, PySAM_dict, pyomo_horizon, 
                       dispatch_time_step):
        """ Initializes the GeneralDispatchParamWrap module
        
        Inputs:
            unitRegistry (pint.registry)   : unique unit Pint unit registry
            SSC_dict (dict)                : dictionary of SSC inputs needed to run modules
            PySAM_dict (dict)              : dictionary of PySAM inputs + file names
            pyomo_horizon (int Quant)      : length of Pyomo simulation segment (hours)
            dispatch_time_step (int Quant) : length of each Pyomo time step (hours)
        """
        
        self.SSC_dict           = SSC_dict
        self.PySAM_dict         = PySAM_dict
        self.pyomo_horizon      = pyomo_horizon
        self.dispatch_time_step = dispatch_time_step 
        
        # setting a standard unit registry
        self.u = unit_registry 
        
        # retrieving the dry bulb temperature values in the SAM directory
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = os.path.join(parent_dir, self.PySAM_dict['solar_resource_rel_parent']) 
        self.Tdry = FileMethods.read_solar_resource_file(self.solar_resource_file, self.u) 
        
        self.set_design()
        

    def set_design(self):
        """ Method to calculate and save design point values of Plant operation
        
        This method extracts values and calculates for design point parameters 
        of our Plant (e.g., nuclear thermal power output, power cycle efficiency,
        inlet and outlet temperatures, etc.). 
        """
        
        u = self.u
        
        # alpha init
        self.alpha = 1.0
        
        # design parameters
        self.p_pb_design  = self.SSC_dict['P_ref'] * u.MW                  # power block design electrical power
        self.eta_design   = self.SSC_dict['design_eff']                    # power block design efficiency
        self.q_pb_design  = (self.p_pb_design / self.eta_design).to('MW')  # power block design thermal rating
        
        # temperature values at design point
        self.T_htf_hot  = (self.SSC_dict['T_htf_hot_des']*u.celsius).to('degK')
        self.T_htf_cold = (self.SSC_dict['T_htf_cold_des']*u.celsius).to('degK')
        
        
    def set_time_indexed_parameters(self, param_dict, df_array, ud_array, current_pyomo_slice):
        """ Method to set time-indexed parameters for Dispatch optimization
        
        This method calculates time parameters for Pyomo Dispatch optimization.
        This includes time-step parameters, which by default are set to 1 hr
        so shouldn't cause problems.
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
            df_array (array)            : array of user defined dispatch factors over simulation time
            ud_array (list of list)     : table of user defined data as nested lists
            current_pyomo_slice (slice) : range of current pyomo horizon (ints representing hours)
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        u = self.u
        
        # time parameters
        pyoH   = self.pyomo_horizon.to('hr').m
        self.T       = int( pyoH )
        self.Delta   = np.array([self.dispatch_time_step.to('hr').m]*self.T)*u.hr
        self.Delta_e = np.cumsum(self.Delta)
        
        # weight parameter
        wt = self.PySAM_dict['weight']
        
        # grab time series data that we have to index
        Tdry  = self.Tdry # dry bulb temperature from solar resource file 
        Price = df_array  # pricing multipliers
        
        # we now curtail any pyomo slice that goes above the full sim time
        length_of_pyoH = current_pyomo_slice.stop - current_pyomo_slice.start
        
        # if we're at the last segment, we won't have 48 hour data for the sim. here is a quick fix
        if length_of_pyoH < pyoH:
            # double-stacking arrays so we can get a slice past full sim time
            Tdry  = np.hstack([Tdry,  Tdry])
            Price = np.hstack([Price, Price])
            
            # restructuring the pyomo slice
            start   = current_pyomo_slice.start
            newstop = start + pyoH
            current_pyomo_slice = slice( start, newstop, 1)
            
        # grabbing relevant dry temperatures
        Tdry   = Tdry[current_pyomo_slice]
        self.P = Price[current_pyomo_slice]*u.USD/u.kWh
        
        # ambient temperature fluctuations and updates
        etamult, wmult = SSCHelperMethods.get_ambient_T_corrections_from_udpc_inputs( self.u, Tdry, ud_array ) # TODO:verify this makes sense
        self.etaamb = etamult * self.SSC_dict['design_eff']
        self.etac   = wmult * self.SSC_dict['ud_f_W_dot_cool_des']/100.

        # power into cycle, into grid, +/-
        self.Qc         = self.Ec / np.ceil(self.SSC_dict['startup_time']*u.hr / np.min(self.Delta)) / np.min(self.Delta) #TODO: make sure Ec is called correctly
        self.Wdotnet    = [1.e10 for j in range(self.T)] *u.kW
        self.W_u_plus   = [(self.Wdotl + self.W_delta_plus*0.5*dt).to('kW').magnitude for dt in self.Delta]*u.kW
        self.W_u_minus  = [(self.Wdotl + self.W_delta_minus*0.5*dt).to('kW').magnitude for dt in self.Delta]*u.kW
        
        # time weights
        n  = len(self.Delta)
        D  = np.zeros(n)
        
        # looping through time
        for t in range(n):
            D[t]        = wt**(self.Delta_e[t]/u.hr)
        self.D = D
        
        #------- Time indexed parameters ---------
        param_dict['T']         = self.T                   #T: time periods
        param_dict['Delta']     = self.Delta.to('hr')      #\Delta_{t}: duration of period t [hr]
        param_dict['Delta_e']   = self.Delta_e.to('hr')    #\Delta_{e,t}: cumulative time elapsed at end of period t [hr]
        param_dict['D']         = self.D                   #D_{t}: Time-weighted discount factor in period $t$ [-]
        param_dict['etaamb']    = self.etaamb              #\eta^{amb}_{t}: Cycle efficiency ambient temperature adjustment factor in period $t$ [-]
        param_dict['etac']      = self.etac                #\eta^{c}_{t}: Normalized condenser parasitic loss in period $t$ [-] 
        param_dict['P']         = self.P.to('USD/kWh')     #P_{t}: Electricity sales price in period $t$ [\$/kWh]
        param_dict['Qc']        = self.Qc.to('kW')         #Q^{c}_{t}: Allowable power per period for cycle start-up in period $t$ [kWt]
        param_dict['Wdotnet']   = self.Wdotnet.to('kW')    #\dot{W}^{net}_{t}: Net grid transmission upper limit in period $t$ [kWe]
        param_dict['W_u_plus']  = self.W_u_plus.to('kW')   #W^{u+}_{t}: Maximum power production when starting generation in period $t$  [kWe]
        param_dict['W_u_minus'] = self.W_u_minus.to('kW')  #W^{u-}_{t}: Maximum power production in period $t$ when stopping generation in period $t+1$  [kWe]
        
        return param_dict

    
    def set_power_cycle_parameters(self, param_dict, ud_array):
        """ Method to set parameters specific to Power Cycle for Dispatch optimization
        
        This method calculates some steady-state parameters specific to the 
        power cycle (PC) for use in Dispatch optimization. These parameters 
        generally have to do with thermal ratings, and upper/lower bounds
        on capacities. 
        
        Inputs:
            param_dict (dict)       : dictionary of Pyomo dispatch parameters
            ud_array (list of list) : table of user defined data as nested lists
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        u = self.u
        # Magic numbers
        Wb_frac = self.PySAM_dict['Wb_frac']             # TODO: get a better estimate- cycle standby parasitic load as frac of cycle capacity
        pc_rampup      = self.PySAM_dict['pc_rampup']      / u.min  # Cycle max ramp up (fraction of capacity per minute)
        pc_rampdown    = self.PySAM_dict['pc_rampdown']    / u.min    
        pc_rampup_vl   = self.PySAM_dict['pc_rampup_vl']   / u.hr   # Cycle ramp up violation limit (fraction of capacity per minute)
        pc_rampdown_vl = self.PySAM_dict['pc_rampdown_vl'] / u.hr 
        
        # fixed parameter calculations
        self.Ec    = self.SSC_dict['startup_frac'] * self.q_pb_design * self.SSC_dict['startup_time']*u.hr      # TODO: this should be an energy, multiplying by min startup time for now
        self.Lc    = (self.SSC_dict['pb_pump_coef']*u.kW/(u.kg/u.s)) * self.dm_pb_design.to('kg/s') / self.q_pb_design.to('kW')
        self.Qb    = self.SSC_dict['q_sby_frac'] * self.q_pb_design
        self.Ql    = self.SSC_dict['cycle_cutoff_frac'] * self.q_pb_design
        self.Qu    = self.SSC_dict['cycle_max_frac'] * self.q_pb_design
        self.Wb    = Wb_frac* self.p_pb_design
        
        #linearized eta calculations
        etap, b = SSCHelperMethods.get_linearized_ud_params( ud_array, self.q_pb_design, self.SSC_dict )
        Wdotl = b + self.Ql*etap
        Wdotu = b + self.Qu*etap
        self.etap  = etap  
        self.Wdotl = Wdotl.to('kW')
        self.Wdotu = Wdotu.to('kW')
        
        # fixed parameter calculations (ctd)
        self.W_delta_plus  = pc_rampup * self.Wdotu 
        self.W_delta_minus = pc_rampdown * self.Wdotu
        self.W_v_plus      = pc_rampup_vl * self.Wdotu
        self.W_v_minus     = pc_rampdown_vl * self.Wdotu
        self.Yu    = self.PySAM_dict['Yu']*u.hr  # TODO: get a better estimate - minimum required power cycle uptime 
        self.Yd    = self.PySAM_dict['Yd']*u.hr  # TODO: get a better estimate - minimum required power cycle downtime 
        
        ### Power Cycle Parameters ###
        param_dict['Ec']      = self.Ec.to('kWh')  #E^c: Required energy expended to start cycle [kWt$\cdot$h]
        param_dict['eta_des'] = self.eta_design    #\eta^{des}: Cycle nominal efficiency [-]
        param_dict['etap']    = self.etap          #\eta^p: Slope of linear approximation of power cycle performance curve [kWe/kWt]
        param_dict['Lc']      = self.Lc.to('')     #L^c: Cycle heat transfer fluid pumping power per unit energy expended [kWe/kWt]
        param_dict['Qb']      = self.Qb.to('kW')   #Q^b: Cycle standby thermal power consumption per period [kWt]
        param_dict['Ql']      = self.Ql.to('kW')   #Q^l: Minimum operational thermal power input to cycle [kWt]
        param_dict['Qu']      = self.Qu.to('kW')   #Q^u: Cycle thermal power capacity [kWt]
        param_dict['Wb']      = self.Wb.to('kW')   #W^b: Power cycle standby operation parasitic load [kWe]
        param_dict['Wdotl']   = self.Wdotl.to('kW')         #\dot{W}^l: Minimum cycle electric power output [kWe]
        param_dict['Wdotu']   = self.Wdotu.to('kW')         #\dot{W}^u: Cycle electric power rated capacity [kWe]
        param_dict['W_delta_plus']  = self.W_delta_plus.to('kW/hr')  #W^{\Delta+}: Power cycle ramp-up designed limit [kWe/h]
        param_dict['W_delta_minus'] = self.W_delta_minus.to('kW/hr') #W^{\Delta-}: Power cycle ramp-down designed limit [kWe/h]
        param_dict['W_v_plus']      = self.W_v_plus.to('kW/hr')      #W^{v+}: Power cycle ramp-up violation limit [kWe/h]
        param_dict['W_v_minus']     = self.W_v_minus.to('kW/hr')     #W^{v-}: Power cycle ramp-down violation limit [kWe/h]
        param_dict['Yu']      = self.Yu.to('hr')   #Y^u: Minimum required power cycle uptime [h]
        param_dict['Yd']      = self.Yd.to('hr')   #Y^d: Minimum required power cycle downtime [h]
        
        return param_dict

    
    def set_fixed_cost_parameters(self, param_dict):
        """ Method to set fixed costs of the Plant
        
        This method calculates some fixed costs for the Plant operations, startup,
        standby, etc. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        u = self.u
        
        # TODO: for now, scaling everything from LORE files
        old_P_ref = self.PySAM_dict['P_ref_baseline'] * u.MW
        P_ratio   = (self.p_pb_design / old_P_ref ).to('')
        
        # TODO: old values from LORE files
        alpha     = self.alpha * u.USD
        C_pc      = self.PySAM_dict['pc_op_cost'] * u.USD/u.kWh        
        C_csu     = self.PySAM_dict['pc_cold_su'] * u.USD
        C_chsp    = self.PySAM_dict['pc_hot_su'] * u.USD
        C_delta_w = self.PySAM_dict['pc_delta_w'] * u.USD/u.kW
        C_v_w     = self.PySAM_dict['pc_delta_w_v'] * u.USD/u.kW
        C_csb     = self.PySAM_dict['pc_sb'] * u.USD/u.kWh

        ### Cost Parameters ###
        param_dict['alpha']       = alpha.to('USD')                    #\alpha: Conversion factor between unitless and monetary values [\$]
        param_dict['Cpc']         = P_ratio * C_pc.to('USD/kWh')       #C^{pc}: Operating cost of power cycle [\$/kWe$\cdot$h]
        param_dict['Ccsu']        = P_ratio * C_csu.to('USD')          #C^{csu}: Penalty for power cycle cold start-up [\$/start]
        param_dict['Cchsp']       = P_ratio * C_chsp.to('USD')         #C^{chsp}: Penalty for power cycle hot start-up [\$/start]
        param_dict['C_delta_w']   = P_ratio * C_delta_w.to('USD/kW')   #C^{\delta_w}: Penalty for change in power cycle  production [\$/$\Delta\text{kWe}$]
        param_dict['C_v_w']       = P_ratio * C_v_w.to('USD/kW')       #C^{v_w}: Penalty for change in power cycle  production tcb{beyond designed limits} [\$/$\Delta\text{kWe}$]
        param_dict['Ccsb']        = P_ratio * C_csb.to('USD/kWh')      #C^{csb}: Operating cost of power cycle standby operation [\$/kWt$\cdot$h]
        
        return param_dict


    def set_initial_state(self, param_dict, updated_dict=None, plant=None, npts=None ):
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
        sschm  = SSCHelperMethods
            
        if updated_dict is None:
            self.current_Plant = copy.deepcopy(self.SSC_dict)
            self.first_run = True
        else:
            self.current_Plant = updated_dict
            self.first_run = False
            
        # important parameters
        e_pb_suinitremain  = self.current_Plant['pc_startup_energy_remain_initial']*u.kWh
        wdot0              = (0 if self.first_run else self.current_Plant['wdot0'])*u.MW 
        y0                 = (self.current_Plant['pc_op_mode_initial'] == 1) 
        ycsb0              = (self.current_Plant['pc_op_mode_initial'] == 2) 
        ycsu0              = (self.current_Plant['pc_op_mode_initial'] == 0 or self.current_Plant['pc_op_mode_initial'] == 4) 
        pc_persist, pc_off = sschm.get_pc_persist_and_off_logs( param_dict, plant, npts ) if plant is not None else [48,48]
        Yu0                = pc_persist if y0       else 0.0
        Yd0                = pc_off     if (not y0) else 0.0
        
        # defining parameters
        self.wdot0 = wdot0.to('kW')  #\dot{w}_0: Initial power cycle electricity generation [kWe] 
        self.y0    = y0              #y_0: 1 if cycle is generating electric power initially, 0 otherwise   
        self.ycsb0 = ycsb0           #y^{csb}_0: 1 if cycle is in standby mode initially, 0 otherwise
        self.ycsu0 = ycsu0           #y^{csu}_0: 1 if cycle is in starting up initially, 0 otherwise
        self.Yu0   = Yu0*u.hr        #Y^u_0: duration that cycle has been generating electric power [h]
        self.Yd0   = Yd0*u.hr        #Y^d_0: duration that cycle has not been generating power (i.e., shut down or in standby mode) [h]
        
        # Initial cycle startup energy accumulated
        tol = 1.e-6
        if np.isnan(e_pb_suinitremain): # SSC seems to report NaN when startup is completed
            self.ucsu0 = self.Ec
        else:   
            self.ucsu0 = max(0.0*u.MWh, self.Ec - e_pb_suinitremain ) 
            if self.ucsu0 > (1.0 - tol)*self.Ec:
                self.ucsu0 = self.Ec
        
        param_dict['ucsu0']  = self.ucsu0.to('kWh')   #u^{csu}_0: Initial cycle start-up energy inventory  [kWt$\cdot$h]
        param_dict['wdot0']  = self.wdot0.to('kW')    #\dot{w}_0: Initial power cycle electricity generation [kW]e
        param_dict['y0']     = self.y0                #y_0: 1 if cycle is generating electric power initially, 0 otherwise
        param_dict['ycsb0']  = self.ycsb0             #y^{csb}_0: 1 if cycle is in standby mode initially, 0 otherwise
        param_dict['ycsu0']  = self.ycsu0             #y^{csu}_0: 1 if cycle is in starting up initially, 0 otherwise
        param_dict['Yu0']    = self.Yu0.to('hr')      #Y^u_0: duration that cycle has been generating electric power [h]
        param_dict['Yd0']    = self.Yd0.to('hr')      #Y^d_0: duration that cycle has not been generating power (i.e., shut down or in standby mode) [h]
        # param_dict['wdot_s_prev']    = 0*u.hr         #\dot{w}^{s,prev}: previous $\dot{w}^s$, or energy sold to grid [kWe]
        # ^ this should be gen[-1] from previous SSC run, 0 if first_run == True
        
        return param_dict
    
    
# =============================================================================
# Dispatch Outputs
# =============================================================================
  
class GeneralDispatchOutputs(object):
    """
    The GeneralDispatchOutputs class is meant to handle outputs from a given,
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
    
