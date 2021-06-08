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
        self.price         = np.asarray( self.mod.TimeOfDeliveryFactors.dispatch_factors_ts )
        
        # mode orders and re-ordering
        op_mode_result, modes_order = np.unique(self.op_mode_1,return_index=True)
        self.op_mode_result = op_mode_result[np.argsort(modes_order)] # re-order modes by first appearance of each
       
        
    def get_array(self, array_str, slicer ):
        # get array from this class object instance
        array = getattr(self, array_str)
        # if the array is in Pint units, grab the magnitude only
        if hasattr(array,'m'):
            array = array.m
        # return the sliced array based on slice input
        return array[slicer]


    def get_slice(self, start_ind, end_ind ):
        # define a slice to use on arrays based on given starting and ending times
        slicer = slice(int(start_ind), int(end_ind), 1)
        return slicer
        

    def plot_on_axis(self, ax, x_array, y_array, label, color=None ):
        # generic plotting given x,y, an axis, label, and optional color for line
        if color is None:
            ax.plot(x_array.m, y_array, linewidth = self.lw, label=label)
        else:
            ax.plot(x_array.m, y_array, color=color, linewidth = self.lw, label=label)


    def bar_plot_on_axis(self, ax, x_array, y_array, dx, label, alpha=0.5, color=None ):
        # generic plotting given x,y, an axis, label, and optional color for line
        if color is None:
            ax.bar(x_array.m, y_array, dx, alpha=0.5, label=label)
        else:
            ax.bar(x_array.m, y_array, dx, color=color, alpha=0.5, label=label)
            
        
    def plot_SSC_generic(self, ax, array_list, label_list, y_label, title_label=None, \
                         plot_all_time=True, start_hr=0, end_hr=48, is_bar_graph=False, return_extra=False ):
        
        # extracting full time array and slice 
        t_plot  = self.t_full.to('d')
        d_slice = self.full_slice
        time_label = 'Time (days)'
        
        # if we're not plotting the full results, slice up the arrays for the time portion we want to plot
        if not plot_all_time:
            d_slice = self.get_slice( start_hr, end_hr )
            t_plot  = self.t_full[d_slice]
            time_label = 'Time (hrs)'
        
        dt = np.diff(t_plot)[0].m
        
        # lambda function to plot arrays to a specific axis
        if is_bar_graph:
            plot_on_axis = lambda axis, array, d_label, color=None:  \
                                  self.bar_plot_on_axis(axis, t_plot, array, dt, d_label, color)
        else:
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
        
        # optional to return the slicer and x-inputs (plotting times)
        if return_extra:
            return ax, d_slice, t_plot
        else:
            return ax
    
    
    def plot_SSC_power_and_energy(self, ax=None, plot_all_time=True, title_label=None, start_hr=0, end_hr=48 ):
        

        # list of array strings 
        power_array_list = [ 'p_cycle' , 'q_dot_rec_in', 'gen']
        
        # list of labels for each array string to extract from Outputs
        power_label_list = [ 'P_cycle (Electric)',
                             'Q_dot to Salt (Thermal)',
                             'Power generated (Electric)' ]
        
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10,5])
            ax = fig.gca()   # this is the power plot
            
        # twin axis to plot energy on opposite y-axis
        ax2 = ax.twinx() # this is the energy plot
        
        # plot Power arrays
        ax = self.plot_SSC_generic(ax, power_array_list, power_label_list, 'Power (MW)', title_label, \
                                            plot_all_time, start_hr, end_hr )
        # plot legend for Power arrays
        ax.legend(loc=self.loc, fontsize=self.fsl)
        
        # custom y limits and ticks to be integers
        ax.set_ylim(0,1100)
        ax.set_yticks([0, 250, 500, 750, 1000])
        
        # custom y limits and ticks to be integers
        ax2.set_ylim(0,1100)
        ax2.set_yticks([0, 250, 500, 750, 1000])

        
        # plot Energy array(s)
        ax2 = self.plot_SSC_generic(ax2, ['e_ch_tes'], ['Salt Charge Level (Thermal)'], 'Energy (MWh)', None, \
                                            plot_all_time, start_hr, end_hr )
        
        # plot legend for Energy arrays and also set line color to default C3 (reddish)
        ax2.get_lines()[0].set_color("C3")
        ax2.legend(loc=self.loc, fontsize=self.fsl)


    def plot_SSC_massflow(self, ax=None, plot_all_time=True, title_label=None, start_hr=0, end_hr=48 ):
        

        # list of array strings 
        power_array_list = [ 'm_dot_pc', 'm_dot_rec'  ]
        
        # list of labels for each array string to extract from Outputs
        power_label_list = [ 'PC HTF mass flow rate', 
                             'Receiver Mass Flow Rate' ]
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10,5])
            ax = fig.gca()   # this is the power plot
            
        # twin axis to plot energy on opposite y-axis
        ax2 = ax.twinx() # this is the defocus plot
        
        # plotting Mass Flow arrays
        ax = self.plot_SSC_generic(ax, power_array_list, power_label_list, 'Mass Flow (kg/s)', title_label, \
                                            plot_all_time, start_hr, end_hr )
        
        # plotting legend for Mass Flow arrays
        ax.legend(loc=self.loc, fontsize=self.fsl)

        # plot Defocus array
        ax2 = self.plot_SSC_generic(ax2, ['defocus'], ['Defocus'], 'Defocus (fraction)', None, \
                                            plot_all_time, start_hr, end_hr )
        
        # plot legend for Energy arrays and also set line color to default C3 (reddish)
        ax2.get_lines()[0].set_color("C3")
        ax2.legend(loc=self.loc, fontsize=self.fsl)
        
        # custom y limits and ticks to be integers
        ax2.set_ylim(0,1.3)
        ax2.set_yticks(np.arange(0,1.1,0.5))


    def plot_SSC_op_modes(self, ax=None, plot_all_time=True, title_label=None, start_hr=0, end_hr=48 ):
        
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10,5])
            ax = fig.gca()   # this is the main plot

        # twin axis to plot energy on opposite y-axis
        ax2 = ax.twinx() # this is the price plot
    
        # create plot for OP Mode line
        ax2 = self.plot_SSC_generic(ax2, ['price'], [' '], 'Price Multiplier ($/kWh)', None, \
                                                       plot_all_time, start_hr, end_hr, True, False )
        
        # custom y limits and ticks to be integers
        ax2.set_ylim(0,2.5)
        ax2.set_yticks(np.arange(0,2.5,0.5))
        
        # create plot for OP Mode line
        ax, d_slice1, t_plot1 = self.plot_SSC_generic(ax, ['op_mode_1'], [' '], 'Operating Mode', title_label, \
                                                       plot_all_time, start_hr, end_hr, False, True )
        
        # plot legend for OP Mode line
        ax.legend(loc=self.loc, fontsize=self.fsl)
        
        # set OP Mode line color to black
        ax.get_lines()[0].set_color("k")
        
        # extract operating modes to designated array
        op_mode_1 = self.get_array('op_mode_1', d_slice1)
        
        # Plotting data points over the OP mode line with different colors and labels
        for op in self.op_mode_result[d_slice1]:
            # getting unique operating modes
            inds = (op_mode_1 == op)
            # individual index getting plotted with unique color and label
            if np.sum(inds):
                ax.plot(t_plot1[inds].m, op_mode_1[inds], 'o', label=self.operating_modes[int(op)])
        
        # plotting legend for OP modes
        ax.legend(loc=self.loc, fontsize=self.fsl)
        
        # custom y limits and ticks to be integers
        ax.set_ylim(0,40)
        ax.set_yticks(np.arange(0,40,5))

        
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

    