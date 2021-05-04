#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 13:40:22 2021

@author: gabrielsoto
"""

import PySAM.GenericSystem as GenericSystem
import PySAM.Grid as Grid
import PySAM.Singleowner as Singleowner
import PySAM.PySSC as pssc
from util.FileMethods import FileMethods
import pint
u = pint.UnitRegistry()
import numpy as np
import copy, os

class GenericSSCModule(object):
    
    def __init__(self, plant_name="generic_system", json_name="100mW_Generic", 
                       is_dispatch=False, dispatch_time_step=1):
        
        # grab names, either default here or from child class
        self.json_name  = json_name
        self.plant_name = plant_name
        
        # read in dictionaries from json script
        PySAM_dict, SSC_dict = FileMethods.read_json( self.json_name )
        
        # storing SSC and Pyomo time horizons, inputs are in unit of hours
        self.ssc_horizon   = PySAM_dict['ssc_horizon'] * u.hr
        self.pyomo_horizon = PySAM_dict['pyomo_horizon'] * u.hr
        self.dispatch_time_step = dispatch_time_step * u.hr
        
        # save SSC_dict for usage later
        self.SSC_dict = SSC_dict
        
        # save csv arrays to class 
        self.store_csv_arrays( PySAM_dict )
        
        # save flag for dispatch
        self.is_dispatch = is_dispatch


    def run_sim(self, run_loop=False):
        """ Method to run single simulation for Generic System
        """
        
        self.run_loop = run_loop
        SSC_dict = self.SSC_dict
        
        #--- create Plant object and execute it
        self.create_Plant( )
        self.simulate_Plant( )
        
        #logging annual energy and capacity factor
        annual_energy = np.sum(self.gen_log)*u.kWh
        if 'P_ref' in SSC_dict.keys():
            ref_energy    = SSC_dict['P_ref']*u.MW * SSC_dict['gross_net_conversion_factor']
            self.capacity_factor = (annual_energy / (ref_energy * 1*u.yr) ).to(' ')
        
        #--- use executed Plant object to create Grid object and execute it
        self.create_Grid( )
        # update gen and annual energy so SystemOutput is consistent, carried over to SO object
        self.Grid.SystemOutput.gen = tuple(self.gen_log)
        self.Grid.SystemOutput.annual_energy = np.sum(annual_energy.magnitude)
        self.Grid.execute( )
        
        #--- use executed Plant object to create SingleOwner object and execute it
        self.create_SO( )
        self.SO.execute( )
        
         
    def store_csv_arrays(self, input_dict):
        """ Method to get data from specified csv files and store in class
        
        Inputs:
            input_dict (dict) : dictionary with csv file names
        """
        
        # saving location of solar resource file for SSC input
        parent_dir = FileMethods.parent_dir
        self.solar_resource_file = os.path.join(parent_dir, input_dict['solar_resource_rel_parent']) #os.path.join
        
        
    def create_Plant(self):
        """ Method to create Plant object for the first time
        """
        
        # create plant data encoding for generic system
        plant_dat = pssc.dict_to_ssc_table( self.SSC_dict, self.plant_name )
        
        # create new Plant object
        self.Plant = GenericSystem.wrap(plant_dat)

    
    def create_Grid(self):
        """ Method to create Grid object for the first time
        """
        
        # create grid data encoding for grid
        grid_dat = pssc.dict_to_ssc_table( self.SSC_dict, "grid" )
        
        # create new Grid object from existing Plant object
        self.Grid = Grid.from_existing( self.Plant )
        
        # import Grid-specific data to Grid object
        self.Grid.assign(Grid.wrap(grid_dat).export())


    def create_SO(self):
        """ Method to create SingleOwner object for the first time
        """
        
        # create singleowner data encoding for singleowner object
        so_dat   = pssc.dict_to_ssc_table( self.SSC_dict, "singleowner" )
        
        # create new Singleowner object from existing Plant object
        self.SO = Singleowner.from_existing( self.Plant )
        
        # import Singleowner-specific data to Singleowner object
        self.SO.assign(Singleowner.wrap(so_dat).export())


    def simulate_Plant(self):
        """ Full simulation of Plant
        """
        
        # helper method for hstack-ing arrays
        self.augment_log = lambda X,Y: np.hstack(  [ X, Y ]  )
        
        # start and end times for full simulation
        time_start = self.SSC_dict['time_start'] * u.s
        time_end   = self.SSC_dict['time_stop'] * u.s
        
        # if running loop -> next 'time stop' is ssc_horizon time
        #            else -> next 'time stop' is sim end time
        time_next  = self.ssc_horizon.to('s') if self.run_loop else copy.deepcopy(time_end)
        
        # setting index for subsequent calls, static index for gen
        self.t_ind = int(time_next.to('hr').magnitude)
        
        # setting up empty log array for gen
        self.gen_log = np.ndarray([0])
        
        # creating parameters for dispatch for first time
        if self.is_dispatch:
            disp_params = self.create_dispatch_params()
        
        # first execution of Plant through SSC
        self.run_Plant_through_SSC( time_start , time_next )
        
        # this loop should only be entered if run_loop == True
        while (time_next < time_end):
            print_time = int(time_next.to('d').magnitude)
            if not print_time % 50: print('   [%s / %s] completed.' % (print_time,time_end.to('d').magnitude))
            
            # update time
            time_start += self.ssc_horizon.to('s')
            time_next  += self.ssc_horizon.to('s')
            
            # run dispatch
            if self.is_dispatch:
                # update dispatch variables from SSC outputs
                disp_vars = self.update_dispatch_vars()
                
                # run dispatch
                outputs = self.run_pyomo()
                
                # update SSC inputs from dispatch outputs
                plant_updt = self.update_Plant_with_dispatch()
            
            # update Plant parameters after previous run
            self.update_Plant( )
            
            # run Plant again
            self.run_Plant_through_SSC( time_start , time_next )
            

    def run_Plant_through_SSC(self, start_hr, end_hr):
        """ Simulation of Plant through SSC for given times
        """
        
        # oops, GenericSystem doesn't have 'SystemControl' subclass...
        if hasattr(self.Plant,'SystemControl'):
            self.Plant.SystemControl.time_start = start_hr.to('s').magnitude
            self.Plant.SystemControl.time_stop = end_hr.to('s').magnitude
        
        self.Plant.execute()
        
        # logging values of gen
        self.gen_log = self.augment_log( self.gen_log, self.Plant.Outputs.gen[0:self.t_ind ] )
        
        
    def reset_all(self):
        """ Reset SSC submodules
        """
        del self.Plant
        del self.Grid
        del self.SO
    
    
    def update_Plant(self):
        """ Replace key Plant inputs with previous outputs in chainlinked simulations
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


    def create_dispatch_params(self):
        
        params = {}

        nhel = 1
        q_rec_design = self.SSC_dict['q_dot_nuclear_des'] * u.MW   # receiver design thermal power
        p_pb_design  = self.SSC_dict['P_ref'] * u.MW               # power block design electrical power
        eta_design   = self.SSC_dict['design_eff']                 # power block efficiency 
        q_pb_design  = p_pb_design / eta_design                    # power block design thermal rating
        
        # TO DOs
        dw_rec_pump  = 0*u.MW             # TODO: Pumping parasitic at design point reciever mass flow rate (MWe)
        tower_piping_ht_loss = 0*u.kW     # TODO: Tower piping heat trace full-load parasitic load (kWe) 
        etap = 0                          # TODO: function needed for this slope
        Wdotl = 0*u.kW                    # TODO: same function as etap
        Wdotu = 0*u.kW                    # TODO: same function as etap
        dm_pb_design = 0*u.kg/u.s         # TODO: get_cycle_design_mass_flow
        Yu = 1*u.hr                       # TODO: minimum required power cycle uptime 
        Yd = 1*u.hr                       # TODO: minimum required power cycle downtime 
        delta_rs = 0                      # TODO: time loop to get this fraction
        D = 0                             # TODO: time loop to get this fraction
        etaamb = 0                        # TODO: function to call ud table to get eta multiplier
        etac = 0                          # TODO: function to call ud table to get wdot multiplier
        
        #------- Time indexed parameters ---------
        params['T']        = int( self.pyomo_horizon.to('hr').magnitude )                       #T: time periods
        params['Delta']    = np.array([self.dispatch_time_step.to('hr').magnitude]*params['T']) #\Delta_{t}: duration of period t
        params['Delta_e']  = np.cumsum(params['Delta'])                                         #\Delta_{e,t}: cumulative time elapsed at end of period t
        
        ### Cost Parameters ###
        params['alpha']  = 1.0
        
        ### CSP Field and Receiver Parameters ###
        params['deltal'] = self.SSC_dict['rec_su_delay']*u.hr      #\delta^l: Minimum time to start the receiver [hr]
        params['Ehs']    = self.SSC_dict['p_start']*u.kWh * nhel   #E^{hs}: Heliostat field startup or shut down parasitic loss [kWe$\cdot$h]
        params['Er']     = (self.SSC_dict['rec_qf_delay'] * q_rec_design).to('kW')   #E^r: Required energy expended to start receiver [kWt$\cdot$h]
        params['Eu']     = (self.SSC_dict['tshours']*u.hr * q_pb_design).to('kWh')   #E^u: Thermal energy storage capacity [kWt$\cdot$h]
        params['Lr']     = (dw_rec_pump / q_rec_design).to('')     #L^r: Receiver pumping power per unit power produced [kWe/kWt]
        params['Qrl']    = (self.SSC_dic['f_rec_min'] * q_rec_design).to('kW')    #Q^{rl}: Minimum operational thermal power delivered by receiver [kWt$\cdot$h]
        params['Qrsb']   = (self.SSC_dic['q_rec_standby_fraction'] * q_rec_design).to('kW')  #Q^{rsb}: Required thermal power for receiver standby [kWt$\cdot$h]
        params['Qrsd']   = (self.SSC_dic['q_rec_shutdown_fraction'] * q_rec_design).to('kW')   #Q^{rsd}: Required thermal power for receiver shut down [kWt$\cdot$h] 
        params['Qru']    = params['Er'] / params['deltal']         #Q^{ru}: Allowable power per period for receiver start-up [kWt$\cdot$h]
        params['Wh']     = self.SSC_dict['p_track']*u.kW           #W^h: Heliostat field tracking parasitic loss [kWe]
        params['Wht']    = tower_piping_ht_loss                    #W^{ht}: Tower piping heat trace parasitic loss [kWe]
        
        ### Power Cycle Parameters ###
        params['Ec']  = (self.SSC_dict['startup_frac'] * q_pb_design).to('kW')  #E^c: Required energy expended to start cycle [kWt$\cdot$h]
        params['eta_des'] = eta_design              #\eta^{des}: Cycle nominal efficiency [-] 
        params['etap']    = etap                    #\eta^p: Slope of linear approximation of power cycle performance curve [kWe/kWt]
        params['Lc']  = ( (self.SSC_dict['pb_pump_coef']*u.kW/u.kg) * dm_pb_design.to('kg') / q_pb_design.to('kW') ).to('')   #L^c: Cycle heat transfer fluid pumping power per unit energy expended [kWe/kWt]
        params['Qb']  = (self.SSC_dict['q_sby_frac'] * q_pb_design).to('kW')         #Q^b: Cycle standby thermal power consumption per period [kWt]
        params['Ql']  = (self.SSC_dict['cycle_cutoff_frac'] * q_pb_design).to('kW')  #Q^l: Minimum operational thermal power input to cycle [kWt]
        params['Qu']  = (self.SSC_dict['cycle_max_frac'] * q_pb_design).to('kW')     #Q^u: Cycle thermal power capacity [kWt]
        params['Wb']  = (self.SSC_dict['Wb_fract']* p_pb_design).to('kW')  #W^b: Power cycle standby operation parasitic load [kWe]
        params['Wdotl'] = Wdotl  #\dot{W}^l: Minimum cycle electric power output [kWe]
        params['Wdotu'] = Wdotu  #\dot{W}^u: Cycle electric power rated capacity [kWe]
        # ramp up/down -> frac/min
        params['W_delta_plus']  = self.SSC_dict['disp_pc_rampup'] * params['Wdotu']      #W^{\Delta+}: Power cycle ramp-up designed limit [kWe/h]
        params['W_delta_minus'] = self.SSC_dict['disp_pc_rampdown'] * params['Wdotu']    #W^{\Delta-}: Power cycle ramp-down designed limit [kWe/h]
        params['W_v_plus']      = self.SSC_dict['disp_pc_rampup_vl'] * params['Wdotu']   #W^{v+}: Power cycle ramp-up violation limit [kWe/h]
        params['W_v_minus']     = self.SSC_dict['disp_pc_rampdown_vl'] * params['Wdotu'] #W^{v-}: Power cycle ramp-down violation limit [kWe/h]
        params['Yu']    = Yu      #Y^u: Minimum required power cycle uptime [h]
        params['Yd']    = Yd      #Y^d: Minimum required power cycle downtime [h]
        
        ### Time series CSP Parameters ###
        params['delta_rs'] = delta_rs #\delta^{rs}_{t}: Estimated fraction of period $t$ required for receiver start-up [-]
        params['D'] = D               #D_{t}: Time-weighted discount factor in period $t$ [-]
        params['etaamb'] = etaamb     #\eta^{amb}_{t}: Cycle efficiency ambient temperature adjustment factor in period $t$ [-]
        params['etac'] = etac         #\eta^{c}_{t}: Normalized condenser parasitic loss in period $t$ [-] 
        params['P']  = self.df_array  #P_{t}: Electricity sales price in period $t$ [\$/kWh]
        params['Qin'] = np.array([q_rec_design.magnitude]*params['T'])*q_rec_design.units   #Q^{in}_{t}: Available thermal power generated by the CSP heliostat field in period $t$ [kWt]
        params['Qc'] = params['Ec'] / np.ceil(self.SSC_dict['startup_time'] / np.min(params['Delta'])) / np.min(params['Delta'])   #Q^{c}_{t}: Allowable power per period for cycle start-up in period $t$ [kWt]
        params['Wdotnet'] = [1.e10 for j in range(params['T'])]  #\dot{W}^{net}_{t}: Net grid transmission upper limit in period $t$ [kWe]
        params['W_u_plus']  = [(params['Wdotl'] + params['W_delta_plus']*0.5*dt) for dt in params['Delta']]   #W^{u+}_{t}: Maximum power production when starting generation in period $t$  [kWe]
        params['W_u_minus'] = [(params['Wdotl'] + params['W_delta_minus']*0.5*dt) for dt in params['Delta']]  #W^{u-}_{t}: Maximum power production in period $t$ when stopping generation in period $t+1$  [kWe]
        
        ### Initial Condition Parameters ###
        
        
        return params
        
        