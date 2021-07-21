#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:44:30 2021

@author: gabrielsoto
"""

import sys
sys.path.append('..')
import PySAM.NuclearTes as NuclearTes
from modules.GenericSSCModule import GenericSSCModule
from dispatch.NuclearDispatch import NuclearDispatch as ND
from dispatch.NuclearDispatch import NuclearDispatchParamWrap as NDP
from dispatch.NuclearDispatch import NuclearDispatchOutputs as NDO
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
import copy

class NuclearTES(GenericSSCModule): 
    """
    The NuclearTES class intializes, updates, and runs SSC simulations through PySAM,
    specifically for the SSC NuclearTES module. 
    """
    
    def __init__(self, plant_name="nuclear_tes", json_name="model1", is_dispatch=False):
        """ Initializes the NuclearTES module
        
        Inputs:
            plant_name (str)         : name of SSC module to run 
            json_name (str)          : name of JSON script with input data for module
            is_dispatch (bool)       : boolean, if True runs Pyomo dispatch optimization
        """
        
        # initialize Generic module, csv data arrays should be saved here
        GenericSSCModule.__init__( self, plant_name, json_name, is_dispatch )
        
        # define specific PySAM module to be called later
        self.PySAM_Module = NuclearTes
        
        # define specific Dispatch module to be called later
        self.Dispatch_Module = ND
        
        self.Dispatch_Outputs = NDO
        
        
    def store_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files and store in class
        
        This method uses the designated PySAM inputs from a JSON script to extract
        csv arrays for use in SSC. The PySAM inputs used here are relative filepaths
        to find the respective csv files. We then either save the filepath as a variable
        or extract the data from the named csv file and save it to as a member attribute
        of this NE2 module class.
        
        Inputs:
            input_dict (dict) : dictionary with csc relative filepaths
        """
        
        # saving location of solar resource file for SSC input using parent class def
        GenericSSCModule.store_csv_arrays(self, input_dict)
        
        # read csv and save data to arrays
        samsim_dir = FileMethods.samsim_dir + '/'
        self.df_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['dispatch_factors_file'])
        self.ud_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['ud_file'])
        self.wl_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['wlim_file'])
        self.hp_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['helio_file'])
        self.gc_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['grid_file'])
        self.em_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['eta_file'])
        self.fm_array = FileMethods.read_csv_through_pandas(samsim_dir + input_dict['flux_file'])
        
        
    def create_Plant(self):
        """ Method to create Plant object for the first time
        
        This method creates a Plant object using built-in PySAM functionalities
        (including some former PySSC structures). Essentially, it creates some
        sort of data structure (pointer?) from SSC inputs found in the SSC_dict
        and the specified SSC module. That data structure is then used to create
        a PySAM module for the specified SSC Plant module (TCSMolten_Salt, etc.). 
        """
        
        # create plant data encoding for generic system
        plant_dat = pssc.dict_to_ssc_table( self.SSC_dict, self.plant_name )
        
        # create new Plant object
        self.Plant = self.PySAM_Module.wrap(plant_dat)
        
        # manually setting data arrays from csv files
        self.Plant.SolarResource.solar_resource_file         = self.solar_resource_file
        self.Plant.TimeOfDeliveryFactors.dispatch_factors_ts = self.df_array
        self.Plant.UserDefinedPowerCycle.ud_ind_od           = self.ud_array
        self.Plant.SystemControl.wlim_series                 = self.wl_array
        self.Plant.HeliostatField.helio_positions            = self.hp_array
        self.Plant.HeliostatField.eta_map                    = self.em_array
        self.Plant.HeliostatField.flux_maps                  = self.fm_array
        self.Plant.SystemControl.dispatch_series = [1.2]*8760


    def create_Grid(self):
        """ Method to create Grid object for the first time
        
        This method creates a Grid object again using built-in PySAM functions.
        The Grid object is created similarly to the Plant object, from SSC inputs
        listed in the SSC_dict. The Grid object, however, is first created
        from the existing Plant object and then the grid-specific input data
        is added to create a wrapper for the SSC Grid module.
        """
        
        # create grid data using parent class
        GenericSSCModule.create_Grid(self)
        
        #set curtailment to be really high
        self.Grid.GridLimits.grid_curtailment = self.gc_array  
        
        
    def run_pyomo(self, params):
        """ Running Pyomo dispatch optimization
        
        ** self.is_dispatch == True
        
        This method strictly runs the Pyomo optimization before execution of an
        SSC segment. It creates a new Dispatch model for the segment, solves it,
        then returns results. Results are stored in a dictionary. 
        
        Inputs:
            params (dict) : dictionary of Pyomo dispatch parameters
        """
        
        # Creation of Dispatch model (could be overloaded)
        dispatch_model = self.Dispatch_Module(params, self.u)
        
        # Solving Dispatch optimization model
        rt_results = dispatch_model.solve_model()
        
        # saving current model to self
        self.current_disp_model = dispatch_model
        
        # retrieving current dispatch counter
        count = str(self.disp_count)
        
        # saving model and results to dictionaries with entries being the current counter value
        self.disp_models[count]  = dispatch_model    
        self.disp_results[count] = rt_results
        
        # increasing dispatch counter value by 1
        self.disp_count += 1



    def create_dispatch_wrapper(self, PySAM_dict):
        """ Creating a wrapper object for calling a class that creates dispatch parameters
        
        ** self.is_dispatch == True 
        (Called in __init__ of NE2 module)
        
        This method creates an object whose class ultimately calculates and creates 
        parameters for Dispatch optimization. The reason this class exists separately
        is that it gets overlaoded based on the PySAM module we are running. Depending on
        the PySAM module, this method calls on a different Dispatch Parameter class that 
        is specific to the module.

        Inputs:
            PySAM_dict (dict)   : dictionary of PySAM inputs from a script in /json directory
        Outputs:
            dispatch_wrap (obj) : wrapper object for the class that creates dispatch parameters
        """
        
        self.DispatchParameterClass = NDP
        
        dispatch_wrap = self.DispatchParameterClass( self.u, self.SSC_dict, PySAM_dict,
                    self.pyomo_horizon, self.dispatch_time_step)
        
        return dispatch_wrap

    
    def create_dispatch_params(self, current_pyomo_slice):
        """ Populating a dictionary with dispatch parameters before optimization
        
        ** self.is_dispatch == True 
        (Called within simulation)
        
        This method is creates the Dispatch Parameter dictionary that will be 
        populated with static inputs from SSC_dict as well as initial conditions
        for Dispatch optimization. The initial conditions are continuously updated
        if simulation is segmented.

        Inputs:
            current_pyomo_slice (slice) : range of current pyomo horizon (ints representing hours)
        Outputs:
            dispatch_wrap (obj) : wrapper object for the class that creates dispatch parameters
        """
        
        # get the object
        DW = self.dispatch_wrap
        
        # run the setters from the GenericSSCModule parent class
        params = GenericSSCModule.create_dispatch_params(self, current_pyomo_slice)
        
        # these are NuclearTES-specific setters
        params = DW.set_nuclear_parameters( params )
        params = DW.set_time_series_nuclear_parameters( params )
        
        # this sets the initial set for the NuclearTES
        params = DW.set_initial_state( params )
        
        return params


    def update_Pyomo_after_SSC(self, params, current_pyomo_slice):
        """ Update Pyomo inputs with SSC outputs from previous segment simulation
        
        ** self.run_loop    == True
        ** self.is_dispatch == True 
                          
        This method uses the SSC end results from the previous simulation segment
        and uses them to update the existing Dispatch parameter dictionary that
        is ultimately sent to Pyomo. Essentially just updates the initial conditions
        of the Dispatch parameter dictionary. 
        
        Inputs:
            params (dict)               : dictionary of Pyomo dispatch parameters
            current_pyomo_slice (slice) : range of current pyomo horizon (ints representing hours)
        Outputs:
            params (dict) : updated dictionary of Pyomo dispatch parameters
        """
        
        updated_SSC_dict = copy.deepcopy(self.SSC_dict)
        
        # saving relevant end-of-sim outputs from the last simulation segment
        updated_SSC_dict['rec_op_mode_initial']              = self.Plant.Outputs.rec_op_mode_final
        updated_SSC_dict['rec_startup_time_remain_init']     = self.Plant.Outputs.rec_startup_time_remain_final
        updated_SSC_dict['rec_startup_energy_remain_init']   = self.Plant.Outputs.rec_startup_energy_remain_final
        updated_SSC_dict['T_tank_cold_init']                 = self.Plant.Outputs.T_tes_cold[self.t_ind-1]
        updated_SSC_dict['T_tank_hot_init']                  = self.Plant.Outputs.T_tes_hot[self.t_ind-1]
        updated_SSC_dict['csp.pt.tes.init_hot_htf_percent']  = self.Plant.Outputs.hot_tank_htf_percent_final
        updated_SSC_dict['pc_op_mode_initial']               = self.Plant.Outputs.pc_op_mode_final
        updated_SSC_dict['pc_startup_time_remain_init']      = self.Plant.Outputs.pc_startup_time_remain_final
        updated_SSC_dict['pc_startup_energy_remain_initial'] = self.Plant.Outputs.pc_startup_energy_remain_final
        
        # these are specific to the initial states
        updated_SSC_dict['wdot0'] = self.Plant.Outputs.P_cycle[self.t_ind-1]
        
        # TODO: removing w_dot_s_prev references in all of Dispatch for now, might need to revisit later
        # updated_SSC_dict['wdot_s_prev'] = 0 #np.array([pe.value(dm.model.wdot_s_prev[t]) for t in dm.model.T])[-1]
        
        DW = self.dispatch_wrap
        
        # updating the initial state and time series Nuclear params
        params = DW.set_time_indexed_parameters( params, self.df_array, self.ud_array, current_pyomo_slice )
        params = DW.set_initial_state( params, updated_SSC_dict, self.Plant, self.t_ind )
        params = DW.set_time_series_nuclear_parameters( params )
        
        return params


    def update_Plant_after_Pyomo(self):
        """ Update SSC Plant inputs with Pyomo optimization outputs from current segment simulation

        ** self.run_loop    == True (can be called outside loop)
        ** self.is_dispatch == True 
        
        This method uses the optimization results from Pyomo and ensures that 
        the next SSC segment uses them throughout the corresponding SSC Horizon.
        SSC normally takes single values for initial conditions (for the first hour
        of the SSC Horizon), but it can also take in an array of values for each
        step in the SSC Horizon. These are called "dispatch_targets". Steps are:
            (1) extract solutions from Pyomo over the Pyomo Horizon, 
            (2) keep the solutions for the shorter SSC Horizon and 
            (3) save these "dispatch target" inputs to the Plant object for the 
                   next SSC simulation segment. 
        """
        
        # number of times in full simulation 
        N_full     = int((self.SSC_dict['time_stop']*self.u.s).to('hr').m)
        
        # the heavy-lifting happens here -> return a dictionary of dispatch target arrays from Pyomo optimization results
        dispatch_targets = self.Dispatch_Outputs.get_dispatch_targets_from_Pyomo( \
                                                self.current_disp_model, self.ssc_horizon, N_full, self.run_loop)
        
        ### Set Dispatch Targets ### 
        # setting dispatch targets to True so that SSC can read in Pyomo inputs
        self.Plant.SystemControl.is_dispatch_targets = True
        
        # extract binary arrays for receiver startup and standby
        self.Plant.SystemControl.is_rec_su_allowed_in = dispatch_targets['is_rec_su_allowed_in']
        self.Plant.SystemControl.is_rec_sb_allowed_in = dispatch_targets['is_rec_sb_allowed_in']
        
        # extract binary arrays for cycle startup and standby
        self.Plant.SystemControl.is_pc_su_allowed_in  = dispatch_targets['is_pc_su_allowed_in']
        self.Plant.SystemControl.is_pc_sb_allowed_in  = dispatch_targets['is_pc_sb_allowed_in']
        
        # extract power arrays for power cycle
        self.Plant.SystemControl.q_pc_target_su_in    = dispatch_targets['q_pc_target_su_in']
        self.Plant.SystemControl.q_pc_target_on_in    = dispatch_targets['q_pc_target_on_in']
        self.Plant.SystemControl.q_pc_max_in          = dispatch_targets['q_pc_max_in']
    