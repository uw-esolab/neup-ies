#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 13:11:44 2021

@author: gabrielsoto
"""

import pint
u = pint.UnitRegistry()
import numpy as np
import matplotlib.pyplot as plt
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)


# BIG difference here is that the Plots class gets initialized
class Plots(object):
    
    def __init__(self, module, lp=16, lps=12, fs=12, lw=2,
                  fsl='x-small', loc='best'):
        
        # full PySAM module
        self.mod = module
        
        # saving full time logs
        self.t_full = np.asarray( self.mod.Outputs.time_hr )*u.hr
        self.full_slice = slice(0, len(self.t_full), 1)
        
        # setting operating modes, kept it way at the bottom because it's ugly
        self.set_operating_modes_list()
        
        # user-defined plotting parameters
        self.lp  = lp   #labelpad
        self.lps = lps  #labelpad short
        self.fs  = fs   #fontsize
        self.lw  = lw   #linewidth
        self.fsl = fsl  #fontsize legend
        self.loc = loc  #location of legend
        
        # alternate legend locations
        self.loc_up = 'upper_right'
        self.loc_lr = 'lower right'   #location of legend
        self.loc_ll = 'lower left'    #location of legend
        self.loc_cr = 'center right'  #location of legend
        self.loc_cl = 'center left'   #location of legend
        
        # saving some outputs for plotting
        self.p_cycle       = np.asarray( self.mod.Outputs.P_cycle )  *u.MW
        self.gen           = (np.asarray( self.mod.Outputs.gen )      *u.kW).to('MW')
        self.q_dot_rec_in  = np.asarray( self.mod.Outputs.q_dot_rec_inc ) *u.MW
        self.m_dot_pc      = np.asarray( self.mod.Outputs.m_dot_pc ) *u.kg/u.s
        self.m_dot_rec     = np.asarray( self.mod.Outputs.m_dot_rec ) *u.kg/u.s
        self.T_pc_in       = np.asarray( self.mod.Outputs.T_pc_in )   *u.degC
        self.T_pc_out      = np.asarray( self.mod.Outputs.T_pc_out )  *u.degC
        self.e_ch_tes      = np.asarray( self.mod.Outputs.e_ch_tes )  *u.MWh
        self.op_mode_1     = np.asarray( self.mod.Outputs.op_mode_1 )
        self.defocus       = np.asarray( self.mod.Outputs.defocus )
        
        # mode orders and re-ordering
        op_mode_result, modes_order = np.unique(self.op_mode_1,return_index=True)
        self.op_mode_result = op_mode_result[np.argsort(modes_order)] # re-order modes by first appearance of each
       
        
    def get_array(self, array_str, slicer ):
        array = getattr(self, array_str)
        if hasattr(array,'m'):
            array = array.m
        return array[slicer]


    def get_slice(self, start_ind, end_ind ):
        slicer = slice(int(start_ind), int(end_ind), 1)
        return slicer
        

    def plot_on_axis(self, ax, x_array, y_array, label, color=None ):
        if color is None:
            ax.plot(x_array, y_array, linewidth = self.lw, label=label)
        else:
            ax.plot(x_array, y_array, color=color, linewidth = self.lw, label=label)
            
        
    def plot_SSC_generic(self, ax, array_list, label_list, y_label, title_label=None, \
                         plot_all_time=True, start_hr=0, end_hr=48, return_extra=False ):
        
        # extracting full time array and slice 
        t_plot  = self.t_full.to('d')
        d_slice = self.full_slice
        time_label = 'Time (days)'
        
        # if we're not plotting the full results, slice up the arrays for the time portion we want to plot
        if not plot_all_time:
            d_slice = self.get_slice( start_hr, end_hr )
            t_plot  = self.t_full[d_slice]
            time_label = 'Time (hrs)'
        
        # lambda function to plot arrays to a specific axis
        plot_on_axis = lambda axis, array, d_label, color=None:  \
                              self.plot_on_axis(axis, t_plot, array, d_label, color)
        
        # lambda function to get arrays from self and slice em
        get_array = lambda array_str: self.get_array( array_str, d_slice )
        
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        # plotting power
        for a,l in zip(array_list, label_list):
            plot_on_axis( ax, get_array(a), l)
        
        # set labels and legend
        ax.set_ylabel(y_label,  labelpad=self.lp, fontsize=self.fs, fontweight='bold')
        ax.set_xlabel(time_label, labelpad=self.lp, fontsize=self.fs, fontweight='bold')
        
        # set title if given
        if title_label is not None:
            ax.set_title(title_label, fontweight='bold')
        
        if return_extra:
            return ax, d_slice, t_plot
        else:
            return ax
    
    
    def plot_SSC_power_and_energy(self, plot_all_time=True, start_hr=0, end_hr=48 ):
        

        # list of array strings 
        power_array_list = [ 'p_cycle' , 'q_dot_rec_in', 'gen']
        
        power_label_list = [ 'P_cycle (Electric)',
                             'Q_dot to Salt (Thermal)',
                             'Power generated (Electric)' ]
        
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        fig = plt.figure(figsize=[10,5])
        ax1 = fig.gca()   # this is the power plot
        ax2 = ax1.twinx() # this is the energy plot
        
        ax1 = self.plot_SSC_generic(ax1, power_array_list, power_label_list, 'Power (MW)', 'SSC_Results', \
                                            plot_all_time, start_hr, end_hr )
        ax1.legend(loc=self.loc, fontsize=self.fsl)
        
        ax2 = self.plot_SSC_generic(ax2, ['e_ch_tes'], ['Salt Charge Level (Thermal)'], 'Energy (MWh)', None, \
                                            plot_all_time, start_hr, end_hr )
            
        ax2.get_lines()[0].set_color("C3")
        ax2.legend(loc=self.loc, fontsize=self.fsl)
        plt.tight_layout()


    def plot_SSC_massflow(self, plot_all_time=True, start_hr=0, end_hr=48 ):
        

        # list of array strings 
        power_array_list = [ 'm_dot_pc', 'm_dot_rec'  ]
        
        power_label_list = [ 'PC HTF mass flow rate', 
                             'Receiver Mass Flow Rate' ]
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        fig = plt.figure(figsize=[10,5])
        ax1 = fig.gca()   # this is the power plot

        ax1 = self.plot_SSC_generic(ax1, power_array_list, power_label_list, 'Mass Flow (kg/s)', 'SSC_Results', \
                                            plot_all_time, start_hr, end_hr )
        ax1.legend(loc=self.loc, fontsize=self.fsl)

        plt.tight_layout()


    def plot_SSC_op_modes(self, plot_all_time=True, start_hr=0, end_hr=48 ):
        
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        fig = plt.figure(figsize=[10,5])
        ax1 = fig.gca()   # this is the main plot
        
        ax1, d_slice1, t_plot1 = self.plot_SSC_generic(ax1, ['op_mode_1'], [' '], 'Operating Mode', 'SSC_Results', \
                                                       plot_all_time, start_hr, end_hr, True )
        ax1.legend(loc=self.loc, fontsize=self.fsl)
        ax1.get_lines()[0].set_color("k")
        op_mode_1 = self.get_array('op_mode_1', d_slice1)
        
        # extra bit: labeling unique operating mode values 
        for op in self.op_mode_result[d_slice1]:
            inds = (op_mode_1 == op)
            if np.sum(inds):
                ax1.plot(t_plot1[inds], op_mode_1[inds], 'o', label=self.operating_modes[int(op)])
        
        ax1.legend(loc=self.loc, fontsize=self.fsl)
        ax1.set_ylim(0,40)
        ax1.set_yticks(np.arange(0,40,5))

        plt.tight_layout()
        
        
    def set_operating_modes_list(self):
        
        # this is taken from ssc/tcs/csp_solver_core.h
        #    in the `C_system_operating_modes` class 
        
        self.operating_modes = [
            'ITER_START',
            'CR_OFF__PC_OFF__TES_OFF__AUX_OFF',
            'CR_SU__PC_OFF__TES_OFF__AUX_OFF',
            'CR_ON__PC_SU__TES_OFF__AUX_OFF',
            'CR_ON__PC_SB__TES_OFF__AUX_OFF',
            'CR_ON__PC_RM_HI__TES_OFF__AUX_OFF',
            'CR_ON__PC_RM_LO__TES_OFF__AUX_OFF',
            'CR_DF__PC_MAX__TES_OFF__AUX_OFF',
            'CR_OFF__PC_SU__TES_DC__AUX_OFF',
            'CR_ON__PC_OFF__TES_CH__AUX_OFF',
            'CR_ON__PC_TARGET__TES_CH__AUX_OFF',
            'CR_ON__PC_TARGET__TES_DC__AUX_OFF',
            'CR_ON__PC_RM_LO__TES_EMPTY__AUX_OFF',
            'CR_DF__PC_OFF__TES_FULL__AUX_OFF',
            'CR_OFF__PC_SB__TES_DC__AUX_OFF',
            'CR_OFF__PC_MIN__TES_EMPTY__AUX_OFF',
            'CR_OFF__PC_RM_LO__TES_EMPTY__AUX_OFF',
            'CR_ON__PC_SB__TES_CH__AUX_OFF',
            'CR_SU__PC_MIN__TES_EMPTY__AUX_OFF',
            'CR_SU__PC_SB__TES_DC__AUX_OFF',
            'CR_ON__PC_SB__TES_DC__AUX_OFF',
            'CR_OFF__PC_TARGET__TES_DC__AUX_OFF',
            'CR_SU__PC_TARGET__TES_DC__AUX_OFF',
            'CR_ON__PC_RM_HI__TES_FULL__AUX_OFF',
            'CR_ON__PC_MIN__TES_EMPTY__AUX_OFF',
            'CR_SU__PC_RM_LO__TES_EMPTY__AUX_OFF',
            'CR_DF__PC_MAX__TES_FULL__AUX_OFF',
            'CR_ON__PC_SB__TES_FULL__AUX_OFF',
            'CR_SU__PC_SU__TES_DC__AUX_OFF',
            'CR_ON__PC_SU__TES_CH__AUX_OFF',
            'CR_DF__PC_SU__TES_FULL__AUX_OFF',
            'CR_DF__PC_SU__TES_OFF__AUX_OFF',
            'CR_TO_COLD__PC_TARGET__TES_DC__AUX_OFF',
            'CR_TO_COLD__PC_RM_LO__TES_EMPTY__AUX_OFF',
            'CR_TO_COLD__PC_SB__TES_DC__AUX_OFF',
            'CR_TO_COLD__PC_MIN__TES_EMPTY__AUX_OFF',
            'CR_TO_COLD__PC_OFF__TES_OFF__AUX_OFF',
            'CR_TO_COLD__PC_SU__TES_DC__AUX_OFF',
            'ITER_END'  ]

        
                
                
            