#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 16:15:50 2021

@author: gabrielsoto
"""

from pylab import rc
import matplotlib.pyplot as plt
import pyomo.environ as pe
import numpy as np
import pint
u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
rc('axes', linewidth=2)
rc('font', weight='bold', size=12)
from util.PostProcessing import Plots

class SolarPlots(Plots):
    """
    The Plots class is a part of the PostProcessing family of classes. It can 
    output results from SSCmodels. This might be a 
    Sisyphus-ian task as some of the parameters are very unique to whatever 
    plot you are trying to make, but will do my best to provide basic structures
    and templates to work off of. 

    Note that the Plots class must be initialized before using. 
    """

    def __init__(self, module, fsl='x-small', loc='best', legend_offset=False,
                 lp=16, lps=12, fs=12, lw=2, x_shrink=0.85, x_legend=12):
        """ Initializes the Plots module

        The instantiation of this class receives a full module object, the module
        being one of the NE2 modules in the /neup-ies/simulations/modules directory.
        It also contains various inputs relating to cosmetic parameters for
        matplotlib plots. 

        Inputs:
            module (object)     : object representing NE2 module class after simulations
            fsl (str)           : fontsize for legend
            loc (str)           : location of legend
            legend_offset(bool) : are we plotting legends off-axis?
            lp (int)            : labelpad for axis labels
            lps (int)           : labelpad for axis labels - short version
            fs (int)            : fontsize for labels, titles, etc.
            lw (int)            : linewidth for plotting
            x_shrink (float)    : (legend_offset==True) amount to shrink axis to make room for legend
        """
        
        Plots.__init__(self, module, fsl=fsl, loc=loc, legend_offset=legend_offset,
                 lp=lp, lps=lps, fs=fs, lw=lw, x_shrink=x_shrink, x_legend=x_legend)
        
        
    def set_ssc_outputs(self, mod_out):
        """ Method to set SSC outputs from module
        
        The method extracts SSC outputs and save them to arrays. The input to this
        method is a subclass to an NE2 class holding all outputs to simulations.          
        
        Inputs:
            mod_out (object)  : object representing Outputs subclass of NE2 module
        """
        # saving some outputs for plotting
        self.p_cycle      = np.asarray(mod_out.P_cycle) * u.MW
        self.gen          = (np.asarray(mod_out.gen) * u.kW).to('MW')
        self.q_dot_rec_in = np.asarray(mod_out.q_dot_rec_inc) * u.MW
        self.q_pb         = np.asarray(mod_out.q_pb) * u.MW
        self.q_dot_pc_su  = np.asarray(mod_out.q_dot_pc_startup) * u.MW
        self.m_dot_pc     = np.asarray(mod_out.m_dot_pc) * u.kg/u.s
        self.m_dot_rec    = np.asarray(mod_out.m_dot_rec) * u.kg/u.s
        self.T_pc_in      = np.asarray(mod_out.T_pc_in) * u.degC
        self.T_pc_out     = np.asarray(mod_out.T_pc_out) * u.degC
        self.T_tes_cold   = np.asarray(mod_out.T_tes_cold) * u.degC
        self.T_tes_hot    = np.asarray(mod_out.T_tes_hot) * u.degC
        self.T_cond_out   = np.asarray(mod_out.T_cond_out) * u.degC
        self.e_ch_tes     = np.asarray(mod_out.e_ch_tes) * u.MWh
        self.op_mode_1    = np.asarray(mod_out.op_mode_1)
        self.defocus      = np.asarray(mod_out.defocus)
        self.price        = np.asarray(self.mod.TimeOfDeliveryFactors.dispatch_factors_ts)

        # setting static inputs
        self.p_pb_design   = self.mod.SystemDesign.P_ref * u.MW              # power block design electrical power
        self.eta_design    = self.mod.SystemDesign.design_eff                # power block design efficiency
        self.q_pb_design  = (self.p_pb_design / self.eta_design).to('MW')    # power block design thermal rating
        self.T_htf_hot    = (self.mod.SystemDesign.T_htf_hot_des*u.celsius).to('degK')   # heat transfer fluid Hot temp
        self.T_htf_cold   = (self.mod.SystemDesign.T_htf_cold_des*u.celsius).to('degK')  # heat transfer fluid Cold temp
        self.e_tes_design = (self.q_pb_design * self.mod.SystemDesign.tshours*u.hr).to('MWh')  # TES storage capacity (kWht)

        # operating modes
        op_mode_result, modes_order = np.unique(self.op_mode_1, return_index=True) # mode orders and re-ordering
        self.op_mode_result = op_mode_result[np.argsort(modes_order)]     # re-order modes by first appearance of each


    def set_operating_modes_list(self):
        """ Method to define list of operating modes

        This method creates a list of operating modes pertaining to the specific SSC module we are
        using. This particular list was taken from /ssc/tcs/csp_solver_core.h, in the 
        `C_system_operating_modes` class.
        """

        self.operating_modes = [
            'ITER_START',
            'CR_OFF__PC_OFF__TES_OFF',
            'CR_SU__PC_OFF__TES_OFF',
            'CR_ON__PC_SU__TES_OFF',
            'CR_ON__PC_SB__TES_OFF',
            'CR_ON__PC_RM_HI__TES_OFF',
            'CR_ON__PC_RM_LO__TES_OFF',
            'CR_DF__PC_MAX__TES_OFF',
            'CR_OFF__PC_SU__TES_DC',
            'CR_ON__PC_OFF__TES_CH',
            'CR_ON__PC_TARGET__TES_CH',
            'CR_ON__PC_TARGET__TES_DC',
            'CR_ON__PC_RM_LO__TES_EMPTY',
            'CR_DF__PC_OFF__TES_FULL',
            'CR_OFF__PC_SB__TES_DC',
            'CR_OFF__PC_MIN__TES_EMPTY',
            'CR_OFF__PC_RM_LO__TES_EMPTY',
            'CR_ON__PC_SB__TES_CH',
            'CR_SU__PC_MIN__TES_EMPTY',
            'CR_SU__PC_SB__TES_DC',
            'CR_ON__PC_SB__TES_DC',
            'CR_OFF__PC_TARGET__TES_DC',
            'CR_SU__PC_TARGET__TES_DC',
            'CR_ON__PC_RM_HI__TES_FULL',
            'CR_ON__PC_MIN__TES_EMPTY',
            'CR_SU__PC_RM_LO__TES_EMPTY',
            'CR_DF__PC_MAX__TES_FULL',
            'CR_ON__PC_SB__TES_FULL',
            'CR_SU__PC_SU__TES_DC',
            'CR_ON__PC_SU__TES_CH',
            'CR_DF__PC_SU__TES_FULL',
            'CR_DF__PC_SU__TES_OFF',
            'CR_TO_COLD__PC_TARGET__TES_DC',
            'CR_TO_COLD__PC_RM_LO__TES_EMPTY',
            'CR_TO_COLD__PC_SB__TES_DC',
            'CR_TO_COLD__PC_MIN__TES_EMPTY',
            'CR_TO_COLD__PC_OFF__TES_OFF',
            'CR_TO_COLD__PC_SU__TES_DC',
            'ITER_END']