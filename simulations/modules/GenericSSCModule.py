#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:40:22 2021

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import PySAM.GenericSystem as GenericSystem
import PySAM.Grid as Grid
import PySAM.Singleowner as Singleowner
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
from util.SSCHelperMethods import SSCHelperMethods
from dispatch.GeneralDispatch import GeneralDispatch as GD
from dispatch.GeneralDispatch import GeneralDispatchParamWrap as GDP
from dispatch.GeneralDispatch import GeneralDispatchOutputs as GDO
from tqdm import tqdm
import numpy as np
import copy
from abc import ABC, abstractmethod

class GenericSSCModule(ABC):
    """
    The GenericSSCModule class works as a way to intialize, update, and
    run SSC simulations through PySAM. 
    
    This class is not meant to be run, but instead work as a parent class
    for all other models including NuclearTES, a future CSP class, etc. 
    
    NOTE: I over-use the word "module" because everything is apparently a module.
    There are three hierarchies of modules at play here:
    - SSC modules   : these are the cmod_<name>.cpp files from the SSC repository written in C++
    - PySAM modules : these are PySAM wrappers for the SSC modules (also Grid, Financial modules, etc.)
    - NE2 modules   : these are Python classes that create PySAM modules (e.g., *THIS* class)
    
    """
    
    @abstractmethod
    def __init__(self, plant_name="abstract", json_name="abstract", 
                       is_dispatch=False, dispatch_time_step=1):
        """ Initializes the GenericSSCModules
        
        Args:
            plant_name (str): 
                name of SSC module to run 
            json_name (str): 
                name of JSON script with input data for module
            is_dispatch (bool): 
                boolean, if True runs Pyomo dispatch optimization
            dispatch_time_step (int): 
                time step for dispatch (hours)
            
        """
        
        # grab names, either default here or from child class
        self.json_name  = json_name
        self.plant_name = plant_name
        
        # create and save a specific unit registry
        self.u = SSCHelperMethods.define_unit_registry()
        u = self.u
        
        # read in dictionaries from json script
        PySAM_dict, SSC_dict, output_keys = FileMethods.read_json( self.json_name )

        # save SSC_dict for usage later
        self.SSC_dict = SSC_dict
        self.PySAM_dict = PySAM_dict
        
        # storing SSC and Pyomo time horizons, inputs are in unit of hours
        self.ssc_horizon   = PySAM_dict['ssc_horizon'] * u.hr
        self.pyomo_horizon = PySAM_dict['pyomo_horizon'] * u.hr
        self.output_keys   = output_keys
        self.dispatch_time_step = dispatch_time_step * u.hr
        
        # save csv arrays to class 
        self.store_csv_arrays( PySAM_dict )
        
        # save flag for dispatch 
        self.is_dispatch = is_dispatch
        
        # initialize dispatch-specific parameters and loggers
        if self.is_dispatch:
            
            # initialize dispatch counter
            self.disp_count = 0
            
            # empty dictionaries to store each dispatch run
            self.disp_models  = {}
            self.disp_results = {}
            
            # initialize dispatch wrap class
            self.dispatch_wrap = self.create_dispatch_wrapper( self.PySAM_dict )


    def run_sim(self, run_loop=False, export=False, filename='temp.csv'):
        """ Method to run single simulation for Generic System
        
        This method handles the creating and execution of Plant, Grid, and Financial objects to
        run through SSC. The Plant has an optional boolean input to allow the running of a full
        simulation in a segmented loop.
        
        Args:
            run_loop (bool): 
                if true, runs simulation in segments. else, runs simulation all at once
            export (bool): 
                if true, exports results to an Excel sheet
            filename (str): 
                name for Excel sheet saved to the /outputs directory
            
        """
        
        self.run_loop = run_loop
        
        #--- create Plant object and execute it
        self.create_Plant( )
        self.simulate_Plant( )
        
        # log final results from looping simulations
        self.log_SSC_arrays(log_final=True)
        
        #--- use executed Plant object to create Grid object and execute it
        self.create_Grid( )
        # update gen and annual energy so SystemOutput is consistent, carried over to SO object
        self.Grid.SystemOutput.gen = tuple(self.gen_log)
        self.Grid.SystemOutput.annual_energy = np.sum(self.annual_energy.magnitude)
        self.Grid.execute( )
        
        #--- use executed Plant object to create SingleOwner object and execute it
        self.create_SO( )
        self.SO.execute( )
        
        if export:
            self.export_results(filename)
         
    @abstractmethod
    def store_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files and store in class
        
        This method uses the designated PySAM inputs from a JSON script to extract
        csv arrays for use in SSC. The PySAM inputs used here are relative filepaths
        to find the respective csv files. We then either save the filepath as a variable
        or extract the data from the named csv file and save it to as a member attribute
        of this NE2 module class.
        
        Args:
            input_dict (dict): 
                dictionary with csc relative filepaths
            
        """
        
        # saving location of solar resource file for SSC input
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = os.path.join(parent_dir, input_dict['solar_resource_rel_parent']) #os.path.join
        
    @abstractmethod    
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
        self.Plant = GenericSystem.wrap(plant_dat)

    @abstractmethod
    def create_Grid(self):
        """ Method to create Grid object for the first time
        
        This method creates a Grid object again using built-in PySAM functions.
        The Grid object is created similarly to the Plant object, from SSC inputs
        listed in the SSC_dict. The Grid object, however, is first created
        from the existing Plant object and then the grid-specific input data
        is added to create a wrapper for the SSC Grid module.
        
        """
        
        # create grid data encoding for grid
        grid_dat = pssc.dict_to_ssc_table( self.SSC_dict, "grid" )
        
        # create new Grid object from existing Plant object
        self.Grid = Grid.from_existing( self.Plant )
        
        # import Grid-specific data to Grid object
        self.Grid.assign(Grid.wrap(grid_dat).export())


    def create_SO(self):
        """ Method to create SingleOwner object for the first time
        
        This method creates a Financial object--in this case, the "SingleOwner"
        SSC financial module-- using PySAM like before. Similarly to the Grid
        object, the SingleOwner object is created first from the existing
        Plant object and then SingleOwner-specific input data is added to create
        a wrapper for the SSC SingleOwner module. 
        
        """
        
        # create singleowner data encoding for singleowner object
        so_dat   = pssc.dict_to_ssc_table( self.SSC_dict, "singleowner" )
        
        # create new Singleowner object from existing Plant object
        self.SO = Singleowner.from_existing( self.Plant )
        
        # import Singleowner-specific data to Singleowner object
        self.SO.assign(Singleowner.wrap(so_dat).export())

    @abstractmethod
    def duplicate_Plant(self, Plant):
        """ Method to create Plant object as a duplicate of existing Plant
        
        This method creates a Plant object from an existing Plant. The new 
        Plant object will have a copy of the original Plant's subclasses
        EXCEPT the Output subclass. The two plant's outputs will NOT be linked.
        
        Note:
            Currently verified in simulations/scripts/sanity_check_scripts
        
        Args:
            Plant (obj): 
                original PySAM Plant module to be copied
        Returns:
            newPlant (obj): 
                duplicate PySAM Plant module, unlinked from original
        """
        
        return Plant
    
    
    def advance_time_segment(self, prev_time_start, prev_time_next):
        """ Method to advance time segment to the next horizon
        

        This method updates the current time segment bookends and slices to
        reflect the next time segment horizon. 
        
        Args:
            prev_time_start (float Quant): 
                starting time of previous time segment
            prev_time_next (float Quant): 
                ending time of previous time segment
        Returns:
            time_start (float Quant): 
                starting time of current SSC horizon
            time_sscH (float Quant): 
                ending time of current SSC horizon
            time_pyoH (float Quant): 
                ending time of current SSC horizon
        """
        
        # advance time segment to next SSC horizon
        time_start = prev_time_start + self.ssc_horizon.to('s')
        time_sscH  = prev_time_next  + self.ssc_horizon.to('s')
        
        # calculate end time for next Pyomo horizon
        time_pyoH = time_start + self.pyomo_horizon.to('s')
        
        # check if the horizon end times exceed the full simulation time
        time_sscH = self.sim_time_end if time_sscH > self.sim_time_end else time_sscH
        time_pyoH = self.sim_time_end if time_pyoH > self.sim_time_end else time_pyoH
        
        # time index for each bookend
        ind_t_start = int( time_start.to('hr').m )
        ind_t_sscH  = int( time_sscH.to('hr').m  )
        ind_t_pyoH  = int( time_pyoH.to('hr').m  )
        
        self.slice_ssc_currentH = slice( ind_t_start, ind_t_sscH, 1)  
        self.slice_pyo_currentH = slice( ind_t_start, ind_t_pyoH, 1) 
        
        return time_start, time_sscH, time_pyoH
        
        
    def simulate_Plant(self):
        """ Method to run full simulation of Plant
        
        This method is a sub-section of the run_sim() method in that it handles the 
        execution of *JUST* the Plant object simulation. However, this constitutes most
        of the SSC computational cost. Namely, this method carries out the setup steps 
        and calls a method to execute the Plant in SSC through the PySAM wrapper. 
        
        Major features:
            (1) "run_loop" (bool) - We can choose to run the simulation all at once, 
                     or run it in segments defined by our "SSC Horizon" defined in the 
                     JSON script.
            (2) "is_dispatch" (bool) - We can also choose to run Dispatch optimization 
                     through a Python package called Pyomo, then use those results to 
                     update our SSC inputs.
        
        The Pyomo optimization is generally conducted over a "Pyomo Horizon" longer than
        the SSC Horizon, but results are only kept for that SSC Horizon to use in the next
        simulation segment. 
        
        """
        
        u = self.u
        
        # start and end times for full simulation
        time_start = self.SSC_dict['time_start'] * u.s
        time_end   = self.SSC_dict['time_stop']  * u.s
        self.sim_time_end = time_end.to('s') 
        
        # determining whether next time segment end is the SSC horizon or full sim time
        time_next  = self.ssc_horizon.to('s') if self.run_loop \
                        else copy.deepcopy(time_end)
        
        # setting index for subsequent calls, it's just a static index for log arrays
        self.t_ind = int(time_next.to('hr').m)
        
        # slice for the first indeces of the SSC Horizon, also initializing current Horizon (updated in loop)
        self.slice_ssc_firstH   = slice( 0, self.t_ind, 1)
        self.slice_pyo_firstH   = slice( 0, int( self.pyomo_horizon.to('hr').m ),  1)
        self.slice_ssc_currentH = slice( int( time_start.to('hr').m ), self.t_ind, 1)  
        self.slice_pyo_currentH = slice( int( time_start.to('hr').m ), int( self.pyomo_horizon.to('hr').m ), 1) 
        
        # setting up empty log array for log arrays
        self.initialize_arrays()
        
        # run dispatch optimization for the first time
        if self.is_dispatch:
            
            # one pre-run of Plant, used to grab inputs to dispatch in some modules
            prePlant = self.duplicate_Plant( self.Plant )
            
            # runs for the full simulation to gather some SSC-specific array calculations
            ssc_run_success, prePlant = self.run_Plant_through_SSC(
                                                prePlant, time_start , time_end 
                                                )
            print("Attempting to run full simulation to gather predictions for Pyomo:  {0}".format(
                            "success!" if ssc_run_success else "failed :("))
            
            # create dispatch parameters for the first time
            disp_params = self.create_dispatch_params( prePlant )
            
            # run pyomo optimization
            self.run_pyomo( disp_params )
            
            # update: Pyomo(t) -> Plant(t) 
            self.Plant = self.update_Plant_after_Pyomo( self.Plant, pre_dispatch_run=False )
            
            del prePlant
        
        # first real execution of Plant through SSC
        ssc_run_success, self.Plant = self.run_Plant_through_SSC( 
                                                self.Plant, time_start , 
                                                time_next 
                                                )
        self.log_SSC_arrays()
        
        # setting up iterable time to cycle through in for loop
        # TODO: this needs to be variable with SSC Horizon
        p_time_next = int(np.round(time_next.to('d').m) )
        p_time_end  = int(np.round(time_end.to('d').m) )
        remaining_sim_time = range(p_time_next, p_time_end) if self.run_loop else range(0)
        
        # this loop should only be entered if run_loop == True
        for t in tqdm(remaining_sim_time):
        
            # advance to the next time segment
            time_start, time_next, time_pyoH = self.advance_time_segment( time_start, time_next )
            
            # update: SSC(t) -> Plant(t+1)
            self.update_Plant_after_SSC( )
            
            # run dispatch optimization
            if self.is_dispatch:
            
                # update: SSC(t) -> Pyomo(t+1)
                disp_params = self.update_Pyomo_after_SSC( self.Plant, disp_params )
                
                # run pyomo optimization again
                self.run_pyomo( disp_params )
            
                # update: Pyomo(t+1) -> Plant(t+1)
                self.Plant = self.update_Plant_after_Pyomo( self.Plant, pre_dispatch_run=False )
            
            # run Plant again
            ssc_run_success, self.Plant = self.run_Plant_through_SSC( \
                                                    self.Plant, time_start , time_next )
            if not ssc_run_success:
                print('\n Plant Simulation ended prematurely.')
                self.log_SSC_arrays( log_final=True )
                break
            
            self.log_SSC_arrays()
        
        if ssc_run_success:
            print('\n Plant Simulation successfully completed.')
        
        
    def run_Plant_through_SSC(self, Plant, start_hr, end_hr):
        """ Simulation of Plant through SSC for given times
        
        This method strictly executes the Plant object in SSC through the PySAM
        wrapper. It updates the start and end times of the simulation in case
        we are running a subsequent segment of the simulation but also handles
        full, un-segmented simulations. If running in segments, it stores outputs
        to member attributes of this NE2 module.
        
        Args:
            Plant (obj): 
                original PySAM Plant module to be executed
            start_hr (float Quant): 
                starting time for next simulation (hours)
            end_hr (float Quant): 
                ending time for next simulation (hours)
        Returns:
            exec_success (bool): 
                True if SSC execution is successful
            Plant (obj): 
                original PySAM Plant module after execution
        """

        # start and end times for full simulation
        Plant.SystemControl.time_start = start_hr.to('s').m
        Plant.SystemControl.time_stop  = end_hr.to('s').m
        
        exec_success = True
        try:
            Plant.execute()
        except Exception as err:
            exec_success = False
            print("\n SSC error: {0}".format(err))
            
        return exec_success, Plant

        
    @abstractmethod
    def run_pyomo(self, params):
        """ Running Pyomo dispatch optimization
        
        Note:
            self.is_dispatch == True
        
        This method strictly runs the Pyomo optimization before execution of an
        SSC segment. It creates a new Dispatch model for the segment, solves it,
        then returns results. Results are stored in a dictionary. 
        
        Args:
            params (dict): 
                dictionary of Pyomo dispatch parameters
            
        """
        
        # Creation of Dispatch model (could be overloaded)
        dispatch_model = GD(params, self.u)
        
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
    
    
    def update_Plant_after_SSC(self):
        """ Update SSC Plant inputs with SSC outputs from previous segment simulation
        
        Note:
            self.run_loop == True
        
        This method uses the SSC end results from the previous simulation segment
        and sets them as the initial conditions for the next SSC segment. As a 
        small note: some outputs are arrays that span the full year, however the
        only relevant parts are the first indeces corresponding to the SSC Horizon.
        All other values are typically 0. 
        
        """
        
        self.Plant.SystemControl.rec_op_mode_initial              = self.Plant.Outputs.rec_op_mode_final
        self.Plant.SystemControl.rec_startup_time_remain_init     = self.Plant.Outputs.rec_startup_time_remain_final
        self.Plant.SystemControl.rec_startup_energy_remain_init   = self.Plant.Outputs.rec_startup_energy_remain_final
        self.Plant.SystemControl.T_tank_cold_init                 = self.Plant.Outputs.T_tes_cold[self.t_ind-1]
        self.Plant.SystemControl.T_tank_hot_init                  = self.Plant.Outputs.T_tes_hot[self.t_ind-1]
        self.Plant.ThermalStorage.csp_pt_tes_init_hot_htf_percent = self.Plant.Outputs.hot_tank_htf_percent_final
        self.Plant.SystemControl.pc_op_mode_initial               = self.Plant.Outputs.pc_op_mode_final
        self.Plant.SystemControl.pc_startup_energy_remain_initial = self.Plant.Outputs.pc_startup_time_remain_final
        self.Plant.SystemControl.pc_startup_time_remain_init      = self.Plant.Outputs.pc_startup_energy_remain_final
      
    @abstractmethod    
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
        return params
        
    @abstractmethod    
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
        # NOTE: will need to overload GDO if this somehow gets overloaded by Nuclear Dispatch Outputs class
        horizon = self.pyomo_horizon if pre_dispatch_run else self.ssc_horizon
        dispatch_targets = GDO.get_dispatch_targets_from_Pyomo(self.current_disp_model, horizon, N_full, self.run_loop)
        
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
        
    @abstractmethod
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
        
        self.DispatchParameterClass = GDP
        
        dispatch_wrap = self.DispatchParameterClass( self.u, self.SSC_dict, PySAM_dict,
                    self.pyomo_horizon, self.dispatch_time_step)
        
        return dispatch_wrap

    @abstractmethod
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
        
        params = {}
        DW = self.dispatch_wrap
        
        # setting parameters for the first time
        params = DW.set_power_cycle_parameters( params, self.ud_array )
        params = DW.set_time_indexed_parameters( params, self.df_array, self.ud_array, self.slice_pyo_firstH )
        params = DW.set_fixed_cost_parameters( params )

        return params


    def initialize_arrays(self):
        """ Initializing empty arrays to log SSC outputs after segment simulations
        
        This method creates empty arrays where SSC outputs will be written to.
        Also creates a list of str names for logged simulation outputs.
        
        """
        
        u = self.u
        
        # start and end times for full simulation
        i_start = (self.SSC_dict['time_start'] * u.s).to('hr').m
        i_end   = (self.SSC_dict['time_stop'] * u.s).to('hr').m
        
        # size of simulation arrays
        N_sim = int( i_end - i_start )
        
        # dictionary of output variable names to log after each segment simulation
        self.Log_Arrays = {
        #    name of NE2 variable || name of SSC module variable
                'time_log':          'time_hr',          # logging time
                'gen_log':           'gen',              # electricity generation log
                'q_thermal_log':     'Q_thermal',        # thermal power from nuclear to HTF 
                'p_cycle_log' :      'P_cycle',          # PC electrical power output (gross)
                'q_dot_rec_inc_log': 'q_dot_rec_inc',    # Nuclear incident thermal power
                'q_pb_log':          'q_pb',             # PC input energy
                'q_dot_pc_su_log' :  'q_dot_pc_startup', # PC startup thermal power
                'm_dot_pc_log' :     'm_dot_pc',         # PC HTF mass flow rate
                'm_dot_rec_log'  :   'm_dot_rec',        # Nuc mass flow rate
                'T_pc_in_log' :      'T_pc_in',          # PC HTF inlet temperature 
                'T_pc_out_log'   :   'T_pc_out',         # PC HTF outlet temperature
                'T_tes_cold_log':    'T_tes_cold',       # TES cold temperature
                'T_tes_hot_log'  :   'T_tes_hot',        # TES hot temperature
                'T_cond_out_log':    'T_cond_out',       # PC condenser water outlet temperature
                'e_ch_tes_log'  :    'e_ch_tes',         # TES charge state
                'op_mode_1_log' :    'op_mode_1',        # Operating Mode
                'defocus_log'   :    'defocus',          # Nuclear "Defocus" fraction
                'eta_log'       :    'eta'               # PC efficiency, gross
            } if self.run_loop \
                 else {'gen_log':    'gen'  # electricity generation log
                      }
        
        # empty array to initalize log arrays
        empty_array = np.zeros(N_sim)
        
        # loop through keys in ^ dictionary, save the KEY name to NE2 module as empty array
        for key in self.Log_Arrays.keys():
            # meta: if we don't grab the copy of empty_array, it'll assign a pointer to the array!!
            setattr( self, key, empty_array.copy() ) 


    def log_SSC_arrays(self, log_final=False):
        """ Creating a wrapper object for calling a class that creates dispatch parameters
        
        This method logs SSC outputs if we choose to run a segmented simulation in a loop. 
        Normally (if log_final == False) it saves outputs to arrays that are members of the 
        GenericSSCModule class; these are previously initialized to 0. (If log_final == True)
        the method gets all the previously logged and filled arrays from the segmented simulation
        and saves it to the Plant object which is a member class of GenericSSCModule. The Plant 
        member class gets a "dummy" subclass made from an empty lambda function but works as a
        class with attributes saved to it. 
        
        Args:
            log_final (bool): 
                if true, save full outputs to Plant. else, log to member arrays of Plant
                
        """
        
        # ==========================================================================
        # wanted to create a quick subclass that where I can extract things during PostProcessing steps
        
        if log_final:
            # don't try this at home...
            self.Plant.PySAM_Outputs = lambda: None
            # quick lambda to make things look pretty
            convert_output = lambda arr: tuple( arr.tolist() )
        
        # ==========================================================================
        # Main Loop -- we're working with pointer magic here, tread carefully
        # ===== eval(Log_Arrays['key']) -> keyword for PySAM_Outputs
        # ================ eval('key')  -> keyword for self
        
        for key in self.Log_Arrays.keys():
            # get current key's value from NE2 module
            self_output  = getattr(self, key )
            
            # if we're still running segmented simulations
            if not log_final:
                # get what we have logged so far
                plant_output = getattr(self.Plant.Outputs,  self.Log_Arrays[key] )
                # grab and save corresponding slices to self (this should be some sort of pointer)
                self_output[self.slice_ssc_currentH] = plant_output[self.slice_ssc_firstH]
            
            # we're done with the full simulation
            else:
                # convert output array to a tuple
                self_output_tuple = convert_output( self_output )
                # save array to the new PySAM_Outputs "subclass"
                setattr( self.Plant.PySAM_Outputs, self.Log_Arrays[key], self_output_tuple )
        
        # ==========================================================================
        # updating capacity factor if logging final OR running non-loop
        
        if log_final or not self.run_loop:
            u = self.u
            # summing up energy for the full simulation and defining reference energy output
            self.annual_energy = np.sum(self.gen_log)*u.kWh
            ref_energy    = self.SSC_dict['P_ref']*u.MW * self.SSC_dict['gross_net_conversion_factor']
            # setting capacity factor to self and (optionally) to PySAM outputs
            self.capacity_factor = (self.annual_energy / (ref_energy * 1*u.yr) ).to(' ')
            
            if log_final:
                setattr( self.Plant.PySAM_Outputs, 'capacity_factor', self.capacity_factor )
                setattr( self.Plant.PySAM_Outputs, 'annual_energy',   self.annual_energy )
                
    
    def export_results(self, filename):
        """ Exports final SSC outputs to a specified Excel file
        
        This method creates an .xlsx file with SSC output data from the full simulation.
        Outputs are specified by keywords in the JSON script supplied by the user. 
        
        Args:
            filename (str): 
                name of xlsx file to save results to within /output directory
        """
        
        # empty lists
        location_list = [] 
        key_list      = []
        value_list    = []
        
        # looping through keywords specified in the JSON script
        for key in self.output_keys['keywords']:
            
            # keyword found in the NE2 class
            if hasattr(self,key):
                location_list.append('PySAM')          # save name of module where keyword found
                value_list.append(getattr(self,key))   # save value of output
            
            # keyword found as an SSC Plant module output
            elif hasattr(self.Plant.Outputs,key):
                location_list.append('Plant')
                value_list.append(getattr(self.Plant.Outputs,key))
            
            # keyword found as an SSC Singleowner module output
            elif hasattr(self.SO.Outputs,key):
                location_list.append('SingleOwner')
                value_list.append(getattr(self.SO.Outputs,key))
            
            # keyword not found, return empty cell
            else:
                location_list.append('Not Found')
                value_list.append(' ')
            
            # append empty list with key name (str)
            key_list.append(key)

        # created nested list with name of module, name of output, and value of output
        dataframe_list = [[x,y,z] for x,y,z in zip(location_list, key_list, value_list)]
        columns=['Object', 'Key', 'Value']
        
        # call the util method to write data to csv/xlsx file
        FileMethods.write_csv(dataframe_list, columns, filename)


    def reset_all(self):
        """ Reset SSC submodules
        
        This method resets all PySAM wrappers, deleting them from this NE2 class.
        Primarily done for unit testing, but could also have use if running 
        simulations in parallel. 
        
        """
        
        def safe_del(x):
            if hasattr(self, x):
                delattr(self, x)
                
        safe_del('Plant')
        safe_del('Grid')
        safe_del('SO')
        
        keys = self.Log_Arrays.keys()
        
        for k in keys:
            safe_del(k)
        
        safe_del('Log_Arrays')