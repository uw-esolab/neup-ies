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
from pyomo.opt import SolverStatus, TerminationCondition
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
        
        Args:
            plant_name (str): 
                name of SSC module to run 
            json_name (str): 
                name of JSON script with input data for module
            is_dispatch (bool): 
                boolean, if True runs Pyomo dispatch optimization
                
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
        
        Args:
            input_dict (dict): 
                dictionary with csv relative filepaths
                
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


    def duplicate_Plant(self, Plant):
        """ Method to create Plant object as a duplicate of existing Plant
        
        This method creates a Plant object from an existing Plant. The new 
        Plant object will have a copy of the original Plant's subclasses
        EXCEPT the Output subclass. The two plant's outputs will NOT be linked.
        
        Note:
            Verified in simulations/scripts/sanity_check_scripts
        
        Args:
            Plant (obj): 
                original PySAM Plant module to be copied
        Returns:
            newPlant (obj): 
                duplicate PySAM Plant module, unlinked from original
                
        """
        # retrieve Plant dictionary
        plant_dict = Plant.export()
        
        # create new Plant object from dictionary
        newPlant = self.PySAM_Module.new()
        newPlant.assign(plant_dict)
        
        return newPlant
        
        
    def run_pyomo(self, params):
        """ Running Pyomo dispatch optimization
        
        Note:
            self.is_dispatch == True
        
        This method strictly runs the Pyomo optimization before execution of an
        SSC segment. It creates a new Dispatch model for the segment, solves it,
        then returns results. Results are stored in a dictionary. Help with failure
        modes was found through here: https://www.pyomo.org/blog/2015/1/8/accessing-solver
        
        Args:
            params (dict): 
                dictionary of Pyomo dispatch parameters
        Returns:
            dispatch_success (bool): 
                if dispatch model was solved successfully, returns True
                
        """
        
        # Creation of Dispatch model (could be overloaded)
        dispatch_model = self.Dispatch_Module(params, self.u)
        
        # Solving Dispatch optimization model
        dispatch_success = True
        try:
            rt_results = dispatch_model.solve_model()
        except Exception as err:
            # usually this gets triggered if cbc solver fails
            dispatch_success = False
            print("Dispatch solver failed with error message: \n{0}".format(err))
        
        # check results to see if saved solution is optimal and feasible (cbc didn't crash)
        if dispatch_success:
            # check if solver status ended normally
            if rt_results.solver.status != SolverStatus.ok:
                dispatch_success = False
                print('Solver solution status was not normal.')
            
            # check if termination condition was not optimal
            if rt_results.solver.termination_condition != TerminationCondition.optimal:
                dispatch_success = False
                print('Solver solution termination condition not optimal: {0}'.format(rt_results.solver.termination_condition) )

        # saving current model to self
        self.current_disp_model = dispatch_model
        
        # retrieving current dispatch counter
        count = str(self.disp_count)
        
        # saving model and results to dictionaries with entries being the current counter value
        self.disp_models[count]  = dispatch_model    
        self.disp_results[count] = rt_results if dispatch_success else {}
        
        # increasing dispatch counter value by 1
        self.disp_count += 1
        
        return dispatch_success


    def create_dispatch_wrapper(self, PySAM_dict):
        """ Creating a wrapper object for calling a class that creates dispatch parameters
        
        Note:
            self.is_dispatch == True 
            (Called in __init__ of NE2 module)
        
        This method creates an object whose class ultimately calculates and creates 
        parameters for Dispatch optimization. The reason this class exists separately
        is that it gets overlaoded based on the PySAM module we are running. Depending on
        the PySAM module, this method calls on a different Dispatch Parameter class that 
        is specific to the module.

        Args:
            PySAM_dict (dict): 
                dictionary of PySAM inputs from a script in /json directory
        Returns:
            dispatch_wrap (obj): 
                wrapper object for the class that creates dispatch parameters
                
        """
        
        self.DispatchParameterClass = NDP
        
        dispatch_wrap = self.DispatchParameterClass( self.u, self.SSC_dict, PySAM_dict,
                    self.pyomo_horizon, self.dispatch_time_step)
        
        return dispatch_wrap

    
    def create_dispatch_params(self, Plant ):
        """ Populating a dictionary with dispatch parameters before optimization
        
        Note:
            self.is_dispatch == True 
            (Called within simulation)
        
        This method is creates the Dispatch Parameter dictionary that will be 
        populated with static inputs from SSC_dict as well as initial conditions
        for Dispatch optimization. The initial conditions are continuously updated
        if simulation is segmented.

        Args:
            Plant (obj): 
                original PySAM Plant module
        Returns:
            dispatch_wrap (obj): 
                wrapper object for the class that creates dispatch parameters
                
        """
        # get the object
        DW = self.dispatch_wrap
        
        # run the setters from the GenericSSCModule parent class
        params = GenericSSCModule.create_dispatch_params(self, Plant )

        # extract array from full run of Plant
        assert hasattr(Plant.Outputs, "Q_thermal"), "Q_thermal was not found in the outputs of Plant."
        self.Q_nuc_guess = Plant.Outputs.Q_thermal
        
        # set up copy of SSC dict
        updated_SSC_dict = copy.deepcopy(self.SSC_dict)
        updated_SSC_dict['Q_thermal'] = self.Q_nuc_guess[self.slice_pyo_firstH]
        
        # these are NuclearTES-specific setters
        params = DW.set_nuclear_parameters( params )
        params = DW.set_time_series_nuclear_parameters( params, updated_SSC_dict )
        
        # this sets the initial set for the NuclearTES
        params = DW.set_initial_state( params )
        
        return params


    def update_Pyomo_after_SSC(self, Plant, params ):
        """ Update Pyomo inputs with SSC outputs from previous segment simulation
        
        Note:
            self.run_loop    == True
            self.is_dispatch == True 
                          
        This method uses the SSC end results from the previous simulation segment
        and uses them to update the existing Dispatch parameter dictionary that
        is ultimately sent to Pyomo. Essentially just updates the initial conditions
        of the Dispatch parameter dictionary. 
        
        Args:
            Plant (obj): 
                original PySAM Plant module
            params (dict): 
                dictionary of Pyomo dispatch parameters
        Returns:
            params (dict): 
                updated dictionary of Pyomo dispatch parameters
                
        """
        
        updated_SSC_dict = copy.deepcopy(self.SSC_dict)
        
        # saving relevant end-of-sim outputs from the last simulation segment
        updated_SSC_dict['rec_op_mode_initial']              = Plant.Outputs.rec_op_mode_final
        updated_SSC_dict['rec_startup_time_remain_init']     = Plant.Outputs.rec_startup_time_remain_final
        updated_SSC_dict['rec_startup_energy_remain_init']   = Plant.Outputs.rec_startup_energy_remain_final
        updated_SSC_dict['T_tank_cold_init']                 = Plant.Outputs.T_tes_cold[self.t_ind-1]
        updated_SSC_dict['T_tank_hot_init']                  = Plant.Outputs.T_tes_hot[self.t_ind-1]
        updated_SSC_dict['csp.pt.tes.init_hot_htf_percent']  = Plant.Outputs.hot_tank_htf_percent_final
        updated_SSC_dict['pc_op_mode_initial']               = Plant.Outputs.pc_op_mode_final
        updated_SSC_dict['pc_startup_time_remain_init']      = Plant.Outputs.pc_startup_time_remain_final
        updated_SSC_dict['pc_startup_energy_remain_initial'] = Plant.Outputs.pc_startup_energy_remain_final
        
        # these are specific to the initial states
        updated_SSC_dict['wdot0'] = Plant.Outputs.P_cycle[self.t_ind-1]

        # extract time series from a previous SSC run
        updated_SSC_dict['Q_thermal'] = self.Q_nuc_guess[self.slice_pyo_currentH]
        
        # TODO: removing w_dot_s_prev references in all of Dispatch for now, might need to revisit later
        # updated_SSC_dict['wdot_s_prev'] = 0 #np.array([pe.value(dm.model.wdot_s_prev[t]) for t in dm.model.T])[-1]
        
        DW = self.dispatch_wrap
        
        # updating the initial state and time series Nuclear params
        params = DW.set_time_indexed_parameters( params, self.df_array, self.ud_array, self.slice_pyo_currentH )
        params = DW.set_initial_state( params, updated_SSC_dict, Plant, self.t_ind )
        params = DW.set_time_series_nuclear_parameters( params, updated_SSC_dict )
        
        return params


    def update_Plant_after_Pyomo(self, Plant, pre_dispatch_run=False):
        """ Update SSC Plant inputs with Pyomo optimization outputs from current segment simulation

        Note:
            self.run_loop    == True (can be called outside loop)
            self.is_dispatch == True 
        
        This method uses the optimization results from Pyomo and ensures that 
        the next SSC segment uses them throughout the corresponding SSC Horizon.
        SSC normally takes single values for initial conditions (for the first hour
        of the SSC Horizon), but it can also take in an array of values for each
        step in the SSC Horizon. These are called "dispatch_targets". Steps are:
        (1) extract solutions from Pyomo over the Pyomo Horizon, 
        (2) keep the solutions for the shorter SSC Horizon and 
        (3) save these "dispatch target" inputs to the Plant object for the 
        next SSC simulation segment. 

        Args:
            Plant (obj): 
                original PySAM Plant module to be updated
            pre_dispatch_run (bool): 
                are we updating the Plant for a pre- or post- dispatch run.
                Recall that we only log post-dispatch Plant runs
        Returns:
            Plant (obj): 
                updated PySAM Plant module
                
        """
        
        # number of times in full simulation 
        N_full     = int((self.SSC_dict['time_stop']*self.u.s).to('hr').m)
        
        # the heavy-lifting happens here -> return a dictionary of dispatch target arrays from Pyomo optimization results
        horizon = self.pyomo_horizon if pre_dispatch_run else self.ssc_horizon
        dispatch_targets = self.Dispatch_Outputs.get_dispatch_targets_from_Pyomo( \
                                                self.current_disp_model, horizon, N_full, self.run_loop)
        
        ### Set Dispatch Targets ### 
        # setting dispatch targets to True so that SSC can read in Pyomo inputs
        Plant.SystemControl.is_dispatch_targets = True
        
        # extract binary arrays for receiver startup and standby
        Plant.SystemControl.is_rec_su_allowed_in = dispatch_targets['is_rec_su_allowed_in']
        Plant.SystemControl.is_rec_sb_allowed_in = dispatch_targets['is_rec_sb_allowed_in']
        
        # extract binary arrays for cycle startup and standby
        Plant.SystemControl.is_pc_su_allowed_in  = dispatch_targets['is_pc_su_allowed_in']
        Plant.SystemControl.is_pc_sb_allowed_in  = dispatch_targets['is_pc_sb_allowed_in']
        
        # extract power arrays for power cycle
        Plant.SystemControl.q_pc_target_su_in    = dispatch_targets['q_pc_target_su_in']
        Plant.SystemControl.q_pc_target_on_in    = dispatch_targets['q_pc_target_on_in']
        Plant.SystemControl.q_pc_max_in          = dispatch_targets['q_pc_max_in']
        
        return Plant
    