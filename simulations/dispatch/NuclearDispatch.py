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
import pint
u = pint.UnitRegistry()

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
    
        # set up costs from parent class
        param_dict = GeneralDispatchParamWrap.set_fixed_cost_parameters( param_dict )
        
        return param_dict
        

    def set_nuclear_parameters(self, param_dict):
        
        # design parameters 
        # TODO: find a way so that these are already set by the time this method is called
        self.q_rec_design = self.SSC_dict['q_dot_nuclear_des'] * u.MW  # receiver design thermal power
        self.p_pb_design  = self.SSC_dict['P_ref'] * u.MW              # power block design electrical power
        self.eta_design   = self.SSC_dict['design_eff']                # power block efficiency
        self.q_pb_design  = self.p_pb_design / self.eta_design         # power block design thermal rating
        
        dw_rec_pump  = 0*u.MW             # TODO: Pumping parasitic at design point reciever mass flow rate (MWe)
        tower_piping_ht_loss = 0*u.kW     # TODO: Tower piping heat trace full-load parasitic load (kWe) 
        
        self.deltal = self.SSC_dict['rec_su_delay']*u.hr
        self.Ehs    = self.SSC_dict['p_start']*u.kWh
        self.Er     = self.SSC_dict['rec_qf_delay'] * self.q_rec_design
        self.Eu     = self.SSC_dict['tshours']*u.hr * self.q_pb_design
        self.Lr     = dw_rec_pump / self.q_rec_design
        self.Qrl    = self.SSC_dict['f_rec_min'] * self.q_rec_design
        self.Qrsb   = self.SSC_dict['q_rec_standby_fraction'] * self.q_rec_design
        self.Qrsd   = self.SSC_dict['q_rec_shutdown_fraction'] * self.q_rec_design
        self.Qru    = self.Er/ self.deltal
        self.Wh     = self.SSC_dict['p_track']*u.kW
        self.Wht    = tower_piping_ht_loss
        
        ### CSP Field and Receiver Parameters ###
        param_dict['deltal'] = self.deltal        #\delta^l: Minimum time to start the receiver [hr]
        param_dict['Ehs']    = self.Ehs.to('kWh') #E^{hs}: Heliostat field startup or shut down parasitic loss [kWe$\cdot$h]
        param_dict['Er']     = self.Er.to('kW')   #E^r: Required energy expended to start receiver [kWt$\cdot$h]
        param_dict['Eu']     = self.Eu.to('kWh')  #E^u: Thermal energy storage capacity [kWt$\cdot$h]
        param_dict['Lr']     = self.Lr.to('')     #L^r: Receiver pumping power per unit power produced [kWe/kWt]
        param_dict['Qrl']    = self.Qrl.to('kW')  #Q^{rl}: Minimum operational thermal power delivered by receiver [kWt$\cdot$h]
        param_dict['Qrsb']   = self.Qrsb.to('kW') #Q^{rsb}: Required thermal power for receiver standby [kWt$\cdot$h]
        param_dict['Qrsd']   = self.Qrsd.to('kW') #Q^{rsd}: Required thermal power for receiver shut down [kWt$\cdot$h] 
        param_dict['Qru']    = self.Qru           #Q^{ru}: Allowable power per period for receiver start-up [kWt$\cdot$h]
        param_dict['Wh']     = self.Wh            #W^h: Heliostat field tracking parasitic loss [kWe]
        param_dict['Wht']    = self.Wht           #W^{ht}: Tower piping heat trace parasitic loss [kWe]
        
        
        return param_dict
    
    
    def set_time_series_nuclear_parameters(self, param_dict):
        pass