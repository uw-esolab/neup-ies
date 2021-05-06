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
from dispatch.GeneralDispatch import GeneralDispatch
from dispatch.GeneralDispatch import GeneralDispatchParamWrap
import pyomo.environ as pe
import numpy as np

class NuclearDispatch(GeneralDispatch):
    
    def __init__(self, params):
        # initialize Generic module, csv data arrays should be saved here
        GeneralDispatch.__init__( self, params )

# =============================================================================
# Dispatch Wrapper
# =============================================================================

class NuclearDispatchParamWrap(GeneralDispatchParamWrap):
    
    def __init__(self, SSC_dict=None, PySAM_dict=None, pyomo_horizon=48, 
                   dispatch_time_step=1):
        
        GeneralDispatchParamWrap.__init__( self, SSC_dict, PySAM_dict, 
                            pyomo_horizon, dispatch_time_step )


    def set_fixed_cost_parameters(self, param_dict):
        
        # grabbing unit registry set up in GeneralDispatch
        u = self.u
    
        # set up costs from parent class
        param_dict = GeneralDispatchParamWrap.set_fixed_cost_parameters( self, param_dict )
        
        C_rec = 0 / u.kWh        
        C_rsu = 0
        C_chsp = 0

        ### Cost Parameters ###
        param_dict['Crec']   = C_rec.to('1/kWh') #C^{rec}: Operating cost of heliostat field and receiver [\$/kWt$\cdot$h]
        param_dict['Crsu']   = C_rsu             #C^{rsu}: Penalty for receiver cold start-up [\$/start]
        param_dict['Crhsp']  = C_chsp            #C^{rhsp}: Penalty for receiver hot start-up [\$/start]

        return param_dict
        

    def set_nuclear_parameters(self, param_dict):
        
        # grabbing unit registry set up in GeneralDispatch
        u = self.u 
        
        time_fix = 1*u.hr                 # TODO: we're missing a time term to fix units
        dw_rec_pump  = 0*u.MW             # TODO: Pumping parasitic at design point reciever mass flow rate (MWe)
        tower_piping_ht_loss = 0*u.kW     # TODO: Tower piping heat trace full-load parasitic load (kWe) 
        q_rec_standby_fraction = 0.05     # TODO: Receiver standby energy consumption (fraction of design point thermal power)
        q_rec_shutdown_fraction = 0.      # TODO: Receiver shutdown energy consumption (fraction of design point thermal power)
        
        self.deltal = self.SSC_dict['rec_su_delay']*u.hr
        self.Ehs    = self.SSC_dict['p_start']*u.kWh
        self.Er     = self.SSC_dict['rec_qf_delay'] * self.q_rec_design * time_fix
        self.Eu     = self.SSC_dict['tshours']*u.hr * self.q_pb_design
        self.Lr     = dw_rec_pump / self.q_rec_design
        self.Qrl    = self.SSC_dict['f_rec_min'] * self.q_rec_design * time_fix
        self.Qrsb   = q_rec_standby_fraction  * self.q_rec_design * time_fix
        self.Qrsd   = q_rec_shutdown_fraction * self.q_rec_design * time_fix
        self.Qru    = self.Er / self.deltal  
        self.Wh     = self.SSC_dict['p_track']*u.kW
        self.Wht    = tower_piping_ht_loss
        
        ### CSP Field and Receiver Parameters ###
        param_dict['deltal'] = self.deltal.to('hr')    #\delta^l: Minimum time to start the receiver [hr]
        param_dict['Ehs']    = self.Ehs.to('kWh')      #E^{hs}: Heliostat field startup or shut down parasitic loss [kWe$\cdot$h]
        param_dict['Er']     = self.Er.to('kWh')       #E^r: Required energy expended to start receiver [kWt$\cdot$h]
        param_dict['Eu']     = self.Eu.to('kWh')       #E^u: Thermal energy storage capacity [kWt$\cdot$h]
        param_dict['Lr']     = self.Lr.to('')          #L^r: Receiver pumping power per unit power produced [kWe/kWt]
        param_dict['Qrl']    = self.Qrl.to('kWh')      #Q^{rl}: Minimum operational thermal power delivered by receiver [kWt$\cdot$h]
        param_dict['Qrsb']   = self.Qrsb.to('kWh')     #Q^{rsb}: Required thermal power for receiver standby [kWt$\cdot$h]
        param_dict['Qrsd']   = self.Qrsd.to('kWh')     #Q^{rsd}: Required thermal power for receiver shut down [kWt$\cdot$h] 
        param_dict['Qru']    = self.Qru.to('kW')      #Q^{ru}: Allowable power per period for receiver start-up [kWt$\cdot$h]
        param_dict['Wh']     = self.Wh.to('kW')        #W^h: Heliostat field tracking parasitic loss [kWe]
        param_dict['Wht']    = self.Wht.to('kW')       #W^{ht}: Tower piping heat trace parasitic loss [kWe]
        
        return param_dict
    
    
    def set_time_series_nuclear_parameters(self, param_dict, df_array):
        
        #MAKE SURE TO CALL THIS METHOD AFTER THE NUCLEAR PARAMETERS 
        u = self.u
        
        delta_rs = 0                      # TODO: time loop to get this fraction
        D = 0                             # TODO: time loop to get this fraction
        etaamb = 0                        # TODO: function to call ud table to get eta multiplier
        etac = 0                          # TODO: function to call ud table to get wdot multiplier
        
        self.delta_rs   = delta_rs
        self.D          = D
        self.etaamb     = etaamb
        self.etac       = etac
        self.P          = df_array
        self.Qin        = np.array([self.q_rec_design.magnitude]*self.T)*self.q_rec_design.units
        self.Qc         = self.Ec / np.ceil(self.SSC_dict['startup_time'] / np.min(self.Delta)) / np.min(self.Delta) #TODO: make sure Ec is called correctly
        self.Wdotnet    = [1.e10 for j in range(self.T)] 
        self.W_u_plus   = [(self.Wdotl + self.W_delta_plus*0.5*dt).to('kW').magnitude for dt in self.Delta]*u.kW
        self.W_u_minus  = [(self.Wdotl + self.W_delta_minus*0.5*dt).to('kW').magnitude for dt in self.Delta]*u.kW
        
        ### Time series CSP Parameters ###
        param_dict['delta_rs']  = self.delta_rs   #\delta^{rs}_{t}: Estimated fraction of period $t$ required for receiver start-up [-]
        param_dict['D']         = self.D          #D_{t}: Time-weighted discount factor in period $t$ [-]
        param_dict['etaamb']    = self.etaamb     #\eta^{amb}_{t}: Cycle efficiency ambient temperature adjustment factor in period $t$ [-]
        param_dict['etac']      = self.etac       #\eta^{c}_{t}: Normalized condenser parasitic loss in period $t$ [-] 
        param_dict['P']         = self.P          #P_{t}: Electricity sales price in period $t$ [\$/kWh]
        param_dict['Qin']       = self.Qin        #Q^{in}_{t}: Available thermal power generated by the CSP heliostat field in period $t$ [kWt]
        param_dict['Qc']        = self.Qc         #Q^{c}_{t}: Allowable power per period for cycle start-up in period $t$ [kWt]
        param_dict['Wdotnet']   = self.Wdotnet    #\dot{W}^{net}_{t}: Net grid transmission upper limit in period $t$ [kWe]
        param_dict['W_u_plus']  = self.W_u_plus   #W^{u+}_{t}: Maximum power production when starting generation in period $t$  [kWe]
        param_dict['W_u_minus'] = self.W_u_minus  #W^{u-}_{t}: Maximum power production in period $t$ when stopping generation in period $t+1$  [kWe]
        
        return param_dict


    def set_initial_state(self, param_dict):
        
        u = self.u
        # can re-use this method by choosing self.SSC_dict if t ==0
        #       or some input dict otherwise
        
        # TES masses, temperatures, specific heat
        m_hot  = self.m_tes_design * (self.SSC_dict['csp.pt.tes.init_hot_htf_percent']/100)        # Available active mass in hot tank
        T_tes_hot_init  = (self.SSC_dict['T_tank_hot_init']*u.celsius).to('degK')
        T_tes_init  = 0.5*(T_tes_hot_init + self.T_htf_cold)
        cp_tes_init = self.get_cp_htf(T_tes_init) 
        
        # important parameters
        e_pb_suinitremain  = self.SSC_dict['pc_startup_energy_remain_initial']*u.kWh
        s_current          = m_hot * cp_tes_init * (T_tes_hot_init - self.T_htf_cold) # TES capacity
        s0                 = min(self.Eu.to('kWh'), s_current.to('kWh')  )
        wdot0              = 0*u.MW
        yr0                = (self.SSC_dict['rec_op_mode_initial'] == 2)
        yrsb0              = False   # TODO: try to use Ty's changes to daotk
        yrsu0              = (self.SSC_dict['rec_op_mode_initial'] == 1)
        y0                 = (self.SSC_dict['pc_op_mode_initial'] == 1) 
        ycsb0              = (self.SSC_dict['pc_op_mode_initial'] == 2) 
        ycsu0              = (self.SSC_dict['pc_op_mode_initial'] == 0 or self.SSC_dict['pc_op_mode_initial'] == 4) 
        disp_pc_persist0   = 0
        disp_pc_off0       = 0
        Yu0                = disp_pc_persist0 if y0 else 0.0
        Yd0                = disp_pc_off0 if (not y0) else 0.0
        Drsu               = 1*u.hr   # Minimum time to start the receiver (hr)
        t_rec_suinitremain = self.SSC_dict['rec_startup_time_remain_init']*u.hr
        e_rec_suinitremain = self.SSC_dict['rec_startup_energy_remain_init']*u.Wh
        rec_accum_time     = max(0.0, Drsu - t_rec_suinitremain )*u.hr
        rec_accum_energy   = max(0.0, self.Er - e_rec_suinitremain )*u.kWh
        # yrsd0             = False   # TODO: try to use Ty's changes to daotk
        # disp_rec_persist0 = 0 
        # drsu0             = disp_rec_persist0 if yrsu0 else 0.0   
        # drsd0             = disp_rec_persist0 if self.SSC_dict['rec_op_mode_initial'] == 0 else 0.0
        
        # defining parameters
        self.s0    = s0              #s_0: Initial TES reserve quantity  [kWt$\cdot$h]
        self.wdot0 = wdot0.to('kW')  #\dot{w}_0: Initial power cycle electricity generation [kWe] 
        self.yr0   = yr0             #y^r_0: 1 if receiver is generating ``usable'' thermal power initially, 0 otherwise 
        self.yrsb0 = yrsb0           #y^{rsb}_0: 1 if receiver is in standby mode initially, 0 otherwise
        self.yrsu0 = yrsu0           #y^{rsu}_0: 1 if receiver is in starting up initially, 0 otherwise
        self.y0    = y0              #y_0: 1 if cycle is generating electric power initially, 0 otherwise   
        self.ycsb0 = ycsb0           #y^{csb}_0: 1 if cycle is in standby mode initially, 0 otherwise
        self.ycsu0 = ycsu0           #y^{csu}_0: 1 if cycle is in starting up initially, 0 otherwise
        self.Yu0   = Yu0             #Y^u_0: duration that cycle has been generating electric power [h]
        self.Yd0   = Yd0             #Y^d_0: duration that cycle has not been generating power (i.e., shut down or in standby mode) [h]
        # self.yrsd0 = yrsd0  # TODO: do we need this? doesn't exist in current GeneralDispatch
        # self.drsu0 = drsu0  # TODO: need this? Duration that receiver has been starting up before the problem horizon (h)
        # self.drsd0 = drsd0  # TODO: defining time in shutdown mode as time "off", will this work in the dispatch model?
        
        # Initial cycle startup energy accumulated
        tol = 1.e-6
        if np.isnan(e_pb_suinitremain): # SSC seems to report NaN when startup is completed
            self.ucsu0 = self.Ec
        else:   
            self.ucsu0 = max(0.0, self.Ec - e_rec_suinitremain ) 
            if self.ucsu0 > (1.0 - tol)*self.Ec:
                self.ucsu0 = self.Ec
        

        # Initial receiver startup energy inventory
        self.ursu0 = min(rec_accum_energy.magnitude, rec_accum_time * self.Qru)  # Note, SS receiver model in ssc assumes full available power is used for startup (even if, time requirement is binding)
        if self.ursu0 > (1.0 - 1.e-6)*self.Er:
            self.ursu0 = self.Er

        # self.ursd0 = 0.0   #TODO: How can we track accumulated shut-down energy (not modeled in ssc)
        
        param_dict['s0']     = self.s0      #s_0: Initial TES reserve quantity  [kWt$\cdot$h]
        param_dict['ucsu0']  = self.ucsu0   #u^{csu}_0: Initial cycle start-up energy inventory  [kWt$\cdot$h]
        param_dict['ursu0']  = self.ursu0   #u^{rsu}_0: Initial receiver start-up energy inventory [kWt$\cdot$h]
        param_dict['wdot0']  = self.wdot0   #\dot{w}_0: Initial power cycle electricity generation [kW]e
        param_dict['yr0']    = self.yr0     #y^r_0: 1 if receiver is generating ``usable'' thermal power initially, 0 otherwise
        param_dict['yrsb0']  = self.yrsb0   #y^{rsb}_0: 1 if receiver is in standby mode initially, 0 otherwise
        param_dict['yrsu0']  = self.yrsu0   #y^{rsu}_0: 1 if receiver is in starting up initially, 0 otherwise
        param_dict['y0']     = self.y0      #y_0: 1 if cycle is generating electric power initially, 0 otherwise
        param_dict['ycsb0']  = self.ycsb0   #y^{csb}_0: 1 if cycle is in standby mode initially, 0 otherwise
        param_dict['ycsu0']  = self.ycsu0   #y^{csu}_0: 1 if cycle is in starting up initially, 0 otherwise
        param_dict['Yu0']    = self.Yu0     #Y^u_0: duration that cycle has been generating electric power [h]
        param_dict['Yd0']    = self.Yd0     #Y^d_0: duration that cycle has not been generating power (i.e., shut down or in standby mode) [h]
        
        return param_dict
