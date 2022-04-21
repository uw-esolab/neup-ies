#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 13:20:58 2022

@author: gabrielsoto
"""

import pyomo.environ as pe
from dispatch.GeneralDispatch import GeneralDispatch
from dispatch.GeneralDispatch import GeneralDispatchParamWrap
from dispatch.IndirectNuclearDispatch import IndirectNuclearDispatch
from dispatch.IndirectNuclearDispatch import IndirectNuclearDispatchParamWrap
from dispatch.SolarDispatch   import SolarDispatch
from dispatch.SolarDispatch   import SolarDispatchParamWrap
import numpy as np
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
import os, copy

class DualIndirectDispatch(SolarDispatch):
    """
    The DualPlantDispatch class is meant to set up and run Dispatch
    optimization as a mixed integer linear program problem using Pyomo,
    specifically for the NuclearMsptTES NE2+SSC module.
    """

    def __init__(self, dual=True, direct=False, **kwargs):
        """ Initializes the IndirectDualDispatch module
        
        The instantiation of this class receives a parameter dictionary from
        the NE2 module (created using the DualPlantDispatchWrapper class). It calls
        on the GeneralDispatch __init__ to create the model. The NuclearDispatch first
        creates an empty Concrete Model from Pyomo, then generates Parameters
        from the parameter dictionary, Variables, Objectives and Constraints.
        
        Inputs:
            params (dict)                : dictionary of Pyomo dispatch parameters
            unitRegistry (pint.registry) : unique unit Pint unit registry
        """
        super().__init__(dual=dual, direct=direct, **kwargs)


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
        super().generate_params(params)
        # generating NuclearDispatch parameters first (PowerCycle, etc.)
        # IndirectNuclearDispatch.generate_params(self, params)
        # SolarDispatch.generate_params(self, params, skip_parent=True)


    def generate_variables(self):
        """ Method to generate parameters within Pyomo DualPlant Model
        
        This method instantiates variables for the Pyomo Concrete Model, with
        domains. Does not need initial guesses here, they are defined in the 
        parameters. We first define continuous and binary variables for the 
        Power Cycle through GeneralDispatch, then declare nuclear variables.
        """
        super().generate_variables()
        # generating NuclearDispatch variables first (PowerCycle, etc.)
        # IndirectNuclearDispatch.generate_variables(self)
        # SolarDispatch.generate_variables(self, skip_parent=True)
    
    
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
                    #obj_cost_nuc_su_hs_sd
                    - (model.Cnsu*model.ynsup[t] + model.Cnhsp*model.ynhsp[t] + model.alpha*model.ynsd[t])
                    #obj_cost_ops
                    - model.Delta[t]*(model.Cpc*model.wdot[t] + model.Ccsb*model.Qb*model.ycsb[t] + model.Crec*model.xr[t] + model.Cnuc*(model.xnp[t] + model.xntes[t]) )
                    for t in model.T ) 
                    )
        
        self.model.OBJ = pe.Objective(rule=objectiveRule, sense = pe.maximize)


    def addPiecewiseLinearEfficiencyConstraints(self):
        """ Method to add efficiency constraints to the Pyomo Solar Model
        
        This method adds constraints pertaining to efficiency constraints defined
        as a piecewise linear approximation. Also referred to as Cycle supply and 
        demand constraints. In the SolarDispatch, we add an extra balance of power
        with respect to energy storage and power produced from the CSP plant. 
        
        TODO: This should be revisited when adding MED!!
        """
        
        super(SolarDispatch, self).addPiecewiseLinearEfficiencyConstraints()
        
        def grid_therm_rule(model, t):
            """ Balance of power flow, i.e. sold vs purchased """
            return (
                    model.wdot_s[t] - model.wdot_p[t] == (1-model.etac[t])*model.wdot[t]
                        - model.Ln*(model.xnp[t] + model.xntes[t] + model.xnsu[t] + model.Qnl*model.ynsb[t])
                		- model.Lr*(model.xr[t] + model.xrsu[t] + model.Qrl*model.yrsb[t])
                		- model.Lc*(model.xnp[t] + model.xtesp[t]) 
                        - model.Wh*model.yr[t] - model.Wb*model.ycsb[t] - model.Wht*(model.yrsb[t]+model.yrsu[t])		#Is Wrsb energy [kWh] or power [kW]?  [az] Wrsb = Wht in the math?
                		- model.Wnht*(model.ynsb[t]+model.ynsu[t])
                        - (model.Ehs/model.Delta[t])*(model.yrsu[t] + model.yrsb[t] + model.yrsd[t])
            )
        
        # call the parent version of this method
        # IndirectNuclearDispatch.addPiecewiseLinearEfficiencyConstraints(self)
        
        # additional constraints
        self.model.del_component( self.model.grid_sun_con )
        self.model.grid_sun_con = pe.Constraint(self.model.T,rule=grid_therm_rule)


    def addTESEnergyBalanceConstraints(self):
        """ Method to add TES constraints to the Pyomo Solar Model
        
        This method adds constraints pertaining to TES energy balance from charging
        with thermal power and discharging to the power cycle. 
        """
        def tes_balance_rule(model, t):
            """ Balance of energy to and from TES """
            if t == 1:
                return model.s[t] - model.s0 == model.Delta[t] * (model.xr[t] + model.xntes[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.xtesp[t] + model.Qnsb*model.ynsb[t] + model.Qrsb*model.yrsb[t]))
            return model.s[t] - model.s[t-1] == model.Delta[t] * (model.xr[t] + model.xntes[t] - (model.Qc[t]*model.ycsu[t] + model.Qb*model.ycsb[t] + model.xtesp[t] + model.Qnsb*model.ynsb[t] + model.Qrsb*model.yrsb[t]))
        def tes_upper_rule(model, t):
            """ Upper bound to TES charge state """
            return model.s[t] <= model.Eu
        def tes_start_up_rule(model, t):
            """ Ensuring sufficient TES charge level to startup NP """
            if t == 1:
                return model.s0 >= model.Delta[t]*model.delta_rs[t]*( (model.Qu + model.Qb)*( -3 + model.yrsu[t] + model.ynsu[t] + model.y0 + model.y[t] + model.ycsb0 + model.ycsb[t] ) + model.xtesp[t] + model.Qb*model.ycsb[t] )
            return model.s[t-1] >= model.Delta[t]*model.delta_rs[t]*( (model.Qu + model.Qb)*( -3 + model.yrsu[t] + model.ynsu[t] + model.y[t-1] + model.y[t] + model.ycsb[t-1] + model.ycsb[t] ) + model.xtesp[t] + model.Qb*model.ycsb[t] )
        def maintain_tes_rule(model):
            """ Final state of TES has to be less than or equal to start """
            return model.s[model.num_periods] <= model.s0
        
        self.model.tes_balance_con = pe.Constraint(self.model.T,rule=tes_balance_rule)
        self.model.tes_upper_con = pe.Constraint(self.model.T,rule=tes_upper_rule)
        self.model.tes_start_up_con = pe.Constraint(self.model.T,rule=tes_start_up_rule)
        self.model.maintain_tes_con = pe.Constraint(rule=maintain_tes_rule)


    def generate_constraints(self):
        """ Method to add ALL constraints to the Pyomo Solar Model
        
        This method calls the previously defined constraint methods to instantiate
        them and add to the existing model. This method first calls the GeneralDispatch
        version to set PowerCycle constraints, then calls nuclear constraint methods
        to add them to the model. 
        """
        
        super().generate_constraints()
        # generating NuclearDispatch constraints first (PowerCycle, etc.)
        # IndirectNuclearDispatch.generate_constraints(self)
        # SolarDispatch.generate_constraints(self, skip_parent=True)

# =============================================================================
# Dispatch Wrapper
# =============================================================================

class DualIndirectDispatchParamWrap(IndirectNuclearDispatchParamWrap):
    """
    The DualPlantDispatchParamWrap class is meant to be the staging area for the 
    creation of Parameters ONLY for the DualPlantDispatch class. It communicates 
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
        
        # here, we invoke the NuclearParamWrap init which calls on the General version
        # within this call, the set_design() method is called which then calls upon both
        #   the nuclear and solar versions
        IndirectNuclearDispatchParamWrap.__init__( self, unit_registry, SSC_dict, PySAM_dict, 
                            pyomo_horizon, dispatch_time_step )


    def set_design(self, skip_parent=False):
        """ Method to calculate and save design point values of Plant operation
        
        This method extracts values and calculates for design point parameters 
        of our Plant (e.g., nuclear thermal power output, power cycle efficiency,
        inlet and outlet temperatures, etc.). 
        """
        
        IndirectNuclearDispatchParamWrap.set_design(self)
        SolarDispatchParamWrap.set_design(self, skip_parent=True, given_des=True)


    def set_fixed_cost_parameters(self, param_dict):
        """ Method to set fixed costs of the Plant
        
        This method calculates some fixed costs for the Plant operations, startup,
        standby, etc. 
        
        Inputs:
            param_dict (dict) : dictionary of Pyomo dispatch parameters
        Outputs:
            param_dict (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        # set up costs from parent class
        param_dict = IndirectNuclearDispatchParamWrap.set_fixed_cost_parameters( self, param_dict )
        param_dict = SolarDispatchParamWrap.set_fixed_cost_parameters(self, param_dict, skip_parent=True) 

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

        # First filling out initial states from GeneralDispatcher
        param_dict = IndirectNuclearDispatchParamWrap.set_initial_state( self, param_dict, updated_dict, plant, npts )
        param_dict = SolarDispatchParamWrap.set_initial_state( self, param_dict, updated_dict, plant, npts, skip_parent=True)
        
        return param_dict


# =============================================================================
# Dispatch Outputs
# =============================================================================
  
class DualIndirectDispatchOutputs(object):
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
        yr   = np.array([pe.value(dm.model.yr[t])   for t in t_pyomo])
        yrsu = np.array([pe.value(dm.model.yrsu[t]) for t in t_pyomo])
        yrsb = np.array([pe.value(dm.model.yrsb[t]) for t in t_pyomo])
        
        #----Cycle Binary Outputs----
        y    = np.array([pe.value(dm.model.y[t])    for t in t_pyomo])
        ycsu = np.array([pe.value(dm.model.ycsu[t]) for t in t_pyomo])
        ycsb = np.array([pe.value(dm.model.ycsb[t]) for t in t_pyomo])
    
        #----Cycle Thermal Power Utilization----
        x = np.array([pe.value(dm.model.xnp[t])   for t in t_pyomo])/1000. # from kWt -> MWt
        x += np.array([pe.value(dm.model.xtesp[t])   for t in t_pyomo])/1000. # from kWt -> MWt
        
        #----Thermal Capacity for Cycle Startup and Operation----
        Qc = np.array([pe.value(dm.model.Qc[t]) for t in t_pyomo])/1000. # from kWt -> MWt
        Qu = dm.model.Qu.value/1000. # from kWt -> MWt
    
        # dispatch target -- nuclear startup/standby binaries
        is_rec_su_allowed_in = [1 if (yr[t] + yrsu[t] + yrsb[t]) > 0.001 else 0 for t in t_horizon]  # Receiver on, startup, or standby
        is_rec_sb_allowed_in = [1 if yrsb[t] > 0.001                     else 0 for t in t_horizon]  # Receiver standby
        
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
            disp_targs['is_rec_su_allowed_in'] = is_rec_su_allowed_in 
            disp_targs['is_rec_sb_allowed_in'] = is_rec_sb_allowed_in
            disp_targs['is_pc_su_allowed_in']  = is_pc_su_allowed_in 
            disp_targs['is_pc_sb_allowed_in']  = is_pc_sb_allowed_in  
            disp_targs['q_pc_target_su_in']    = q_pc_target_su_in  
            disp_targs['q_pc_target_on_in']    = q_pc_target_on_in
            disp_targs['q_pc_max_in']          = q_pc_max_in
        # if we're running full simulation all at once, need arrays to match size of full sim
        # TODO: is this a feature we want in the long term? Or just for debugging the first Pyomo call?
        else:
            disp_targs['is_rec_su_allowed_in'] = np.hstack( [is_rec_su_allowed_in , empty_array] ).tolist()
            disp_targs['is_rec_sb_allowed_in'] = np.hstack( [is_rec_sb_allowed_in , empty_array] ).tolist()
            disp_targs['is_pc_su_allowed_in']  = np.hstack( [is_pc_su_allowed_in  , empty_array] ).tolist()
            disp_targs['is_pc_sb_allowed_in']  = np.hstack( [is_pc_sb_allowed_in  , empty_array] ).tolist()
            disp_targs['q_pc_target_su_in']    = np.hstack( [q_pc_target_su_in    , empty_array] ).tolist()
            disp_targs['q_pc_target_on_in']    = np.hstack( [q_pc_target_on_in    , empty_array] ).tolist()
            disp_targs['q_pc_max_in']          = np.hstack( [q_pc_max_in          , empty_array] ).tolist()
            
        return disp_targs