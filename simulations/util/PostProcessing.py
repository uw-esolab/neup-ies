#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 13:11:44 2021

@author: gabrielsoto
"""

import pint
u = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
import numpy as np
import pyomo.environ as pe
import matplotlib.pyplot as plt
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)


# BIG difference here is that the Plots class gets initialized
class Plots(object):
    
    def __init__(self, module, lp=16, lps=12, fs=12, lw=2,
                  fsl='x-small', loc='best'):
        
        # full PySAM module
        self.mod = module.Plant
        
        # define an Output object to extract information from
        Outputs = self.mod.PySAM_Outputs if module.run_loop else self.mod.Outputs
        
        # saving full time logs
        self.t_full = np.asarray( Outputs.time_hr )*u.hr
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
        self.loc_ur = 'upper right'
        self.loc_ul = 'upper left'
        self.loc_lr = 'lower right'   #location of legend
        self.loc_ll = 'lower left'    #location of legend
        self.loc_cr = 'center right'  #location of legend
        self.loc_cl = 'center left'   #location of legend
        
        # saving some outputs for plotting
        self.p_cycle       = np.asarray( Outputs.P_cycle )  *u.MW
        self.gen           = (np.asarray( Outputs.gen )      *u.kW).to('MW')
        self.q_dot_rec_in  = np.asarray( Outputs.q_dot_rec_inc ) *u.MW
        self.q_pb          = np.asarray( Outputs.q_pb ) *u.MW
        self.q_dot_pc_su   = np.asarray( Outputs.q_dot_pc_startup ) *u.MW
        self.m_dot_pc      = np.asarray( Outputs.m_dot_pc )   *u.kg/u.s
        self.m_dot_rec     = np.asarray( Outputs.m_dot_rec )  *u.kg/u.s
        self.T_pc_in       = np.asarray( Outputs.T_pc_in )    *u.degC  
        self.T_pc_out      = np.asarray( Outputs.T_pc_out )   *u.degC
        self.T_tes_cold    = np.asarray( Outputs.T_tes_cold ) *u.degC
        self.T_tes_hot     = np.asarray( Outputs.T_tes_hot )  *u.degC
        self.T_cond_out    = np.asarray( Outputs.T_cond_out )  *u.degC
        self.e_ch_tes      = np.asarray( Outputs.e_ch_tes )   *u.MWh
        self.op_mode_1     = np.asarray( Outputs.op_mode_1 )
        self.defocus       = np.asarray( Outputs.defocus )
        self.price         = np.asarray( self.mod.TimeOfDeliveryFactors.dispatch_factors_ts )
        
        # static inputs
        self.q_rec_design = self.mod.SystemDesign.q_dot_nuclear_des * u.MW      # receiver design thermal power
        self.p_pb_design  = self.mod.SystemDesign.P_ref * u.MW                  # power block design electrical power
        self.eta_design   = self.mod.SystemDesign.design_eff                   # power block design efficiency
        self.q_pb_design  = (self.p_pb_design / self.eta_design).to('MW')  # power block design thermal rating
        self.T_htf_hot  = (self.mod.SystemDesign.T_htf_hot_des*u.celsius).to('degK')
        self.T_htf_cold = (self.mod.SystemDesign.T_htf_cold_des*u.celsius).to('degK')
        self.e_tes_design = (self.q_pb_design * self.mod.SystemDesign.tshours*u.hr).to('MWh') # TES storage capacity (kWht)
        
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
        power_array_list = [ 'p_cycle' , 'q_dot_rec_in', 'gen', 'q_dot_pc_su']
        
        # list of labels for each array string to extract from Outputs
        power_label_list = [ 'P_cycle (Electric)',
                             'Q_dot to Salt (Thermal)',
                             'Power generated (Electric)',
                             'PC startup thermal power (Thermal)']
        
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
        
        # custom y limits and ticks to be integers
        ax2.set_ylim(0,5000)
        ax2.set_yticks(np.linspace(0,5000,5) )
        # ax2.set_ylim(-0.05*self.e_tes_design.m,self.e_tes_design.m)
        # ax2.set_yticks( np.round( np.linspace(0,self.e_tes_design.m, 5), -4 ) )

        # plot Energy array(s)
        ax2 = self.plot_SSC_generic(ax2, ['e_ch_tes'], ['Salt Charge Level (Thermal)'], 'Energy (MWh)', None, \
                                            plot_all_time, start_hr, end_hr )
        
        # plot legend for Power arrays
        ax.legend(loc=self.loc_ul, fontsize=self.fsl)
        
        # custom y limits and ticks to be integers
        ax.set_ylim(-100,1100)
        ax.set_yticks([0, 250, 500, 750, 1000])
        
        # plot legend for Energy arrays and also set line color to default C3 (reddish)
        ax2.get_lines()[0].set_color("C3")
        ax2.legend(loc=self.loc_ur, fontsize=self.fsl)


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


    def plot_SSC_temperatures(self, ax=None, plot_all_time=True, title_label=None, start_hr=0, end_hr=48 ):
        
        # list of array strings 
        power_array_list = [ 'T_pc_in', 'T_pc_out', 'T_tes_cold', 'T_tes_hot'  ]
        
        # list of labels for each array string to extract from Outputs
        power_label_list = [ 'PC HTF (hot) inlet temperature', 
                             'PC HTF (cold) outlet temperature',
                             'TES cold temperature',
                             'TES hot temperature']
        
        #-------------------------#
        #---- Creating Figure ----#
        #-------------------------#
        
        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10,5])
            ax = fig.gca()   # this is the power plot
            
        # plotting Mass Flow arrays
        ax = self.plot_SSC_generic(ax, power_array_list, power_label_list, 'Temperature (C)', title_label, \
                                            plot_all_time, start_hr, end_hr )
        
        # plotting legend for Mass Flow arrays
        ax.legend(loc=self.loc, fontsize=self.fsl)

        

    def plot_pyomo(self, dm):
        
        lw = self.lw
        lp = self.lp
        lps = self.lps
        fs = self.fs
        fsl = self.fsl
        loc = 'best'
        loc_cr = self.loc_cr
        
        #___Time Array
        t_array = np.array([pe.value(dm.model.Delta_e[t]) for t in dm.model.T])
        
        #__Price Array
        p_array = np.array([pe.value(dm.model.P[t]) for t in dm.model.T])
        d_array = np.array([pe.value(dm.model.D[t]) for t in dm.model.T])
        
        #__Costs
        Ccsu      = pe.value(dm.model.Ccsu)
        Cchsp     = pe.value(dm.model.Cchsp)
        C_delta_w = pe.value(dm.model.C_delta_w)
        C_v_w     = pe.value(dm.model.C_v_w)
        Crsu      = pe.value(dm.model.Crsu)
        Crhsp     = pe.value(dm.model.Cchsp)
        Cpc       = pe.value(dm.model.Cpc)
        Ccsb      = pe.value(dm.model.Ccsb)
        Crec      = pe.value(dm.model.Crec)
        Qb        = pe.value(dm.model.Qb)
        
        #___Energy Arrays
        s_array = ( np.array([pe.value(dm.model.s[t]) for t in dm.model.T])*u.kWh ).to('MWh')
        ucsu_array = ( np.array([pe.value(dm.model.ucsu[t]) for t in dm.model.T])*u.kWh ).to('MWh')
        ursu_array = ( np.array([pe.value(dm.model.ursu[t]) for t in dm.model.T])*u.kWh ).to('MWh')
        
        #___Power Arrays
        wdot_array = ( np.array([pe.value(dm.model.wdot[t]) for t in dm.model.T])*u.kW ).to('MW')
        wdot_delta_plus_array = ( np.array([pe.value(dm.model.wdot_delta_plus[t]) for t in dm.model.T])*u.kW ).to('MW')
        wdot_delta_minus_array = ( np.array([pe.value(dm.model.wdot_delta_minus[t]) for t in dm.model.T])*u.kW ).to('MW')
        wdot_v_plus_array = ( np.array([pe.value(dm.model.wdot_v_plus[t]) for t in dm.model.T])*u.kW ).to('MW')
        wdot_v_minus_array = ( np.array([pe.value(dm.model.wdot_v_minus[t]) for t in dm.model.T])*u.kW ).to('MW')
        wdot_s_array = ( np.array([pe.value(dm.model.wdot_s[t]) for t in dm.model.T])*u.kW ).to('MW')
        wdot_p_array = ( np.array([pe.value(dm.model.wdot_p[t]) for t in dm.model.T])*u.kW ).to('MW')
        x_array = ( np.array([pe.value(dm.model.x[t]) for t in dm.model.T])*u.kW ).to('MW')
        xr_array = ( np.array([pe.value(dm.model.xr[t]) for t in dm.model.T])*u.kW ).to('MW')
        xrsu_array = ( np.array([pe.value(dm.model.xrsu[t]) for t in dm.model.T])*u.kW ).to('MW')
        # wdot_s_prev_delta_plus_array = ( np.array([pe.value(dm.model.wdot_s_prev_delta_plus[t]) for t in dm.model.T])*u.kW ).to('MW')
        # wdot_s_prev_delta_minus_array = ( np.array([pe.value(dm.model.wdot_s_prev_delta_minus[t]) for t in dm.model.T])*u.kW ).to('MW')
        
        #___Binary Arrays
        yr_array = np.array([pe.value(dm.model.yr[t]) for t in dm.model.T])
        yrhsp_array = np.array([pe.value(dm.model.yrhsp[t]) for t in dm.model.T])
        yrsb_array = np.array([pe.value(dm.model.yrsb[t]) for t in dm.model.T])
        yrsd_array = np.array([pe.value(dm.model.yrsd[t]) for t in dm.model.T])
        yrsu_array = np.array([pe.value(dm.model.yrsu[t]) for t in dm.model.T])
        yrsup_array = np.array([pe.value(dm.model.yrsup[t]) for t in dm.model.T])
        y_array = np.array([pe.value(dm.model.y[t]) for t in dm.model.T])+2
        ychsp_array = np.array([pe.value(dm.model.ychsp[t]) for t in dm.model.T])+2
        ycsb_array = np.array([pe.value(dm.model.ycsb[t]) for t in dm.model.T])+2
        ycsd_array = np.array([pe.value(dm.model.ycsd[t]) for t in dm.model.T])+2
        ycsu_array = np.array([pe.value(dm.model.ycsu[t]) for t in dm.model.T])+2
        ycsup_array = np.array([pe.value(dm.model.ycsup[t]) for t in dm.model.T])+2
        ycgb_array = np.array([pe.value(dm.model.ycgb[t]) for t in dm.model.T])+2
        ycge_array = np.array([pe.value(dm.model.ycge[t]) for t in dm.model.T])+2
        
        #___Marking 24 hour lines
        time_24hr   = np.ones([len(t_array)])*24
        energy_24hr = np.linspace( np.min( s_array.m ),
                                   np.max( s_array.m, )*1.1, len(t_array) )
        power_24hr  = np.linspace( np.min( [wdot_array.m, x_array.m] ),
                                   np.max( [wdot_array.m, x_array.m] )*1.3, len(t_array) )
        binary_24hr = np.linspace(0,3.1,len(t_array))
        binary_div  = 1.5*np.ones(len(t_array))
        
        # ======= Plotting ========== #
        fig1 = plt.figure(figsize=[12,10])
        ax1 = fig1.add_subplot(411)
        ax2 = fig1.add_subplot(412)
        ax3 = fig1.add_subplot(413)
        ax4 = fig1.add_subplot(414)
        ax5 = ax2.twinx()
        
        #___Shrinking x-axis to allow room for legends off-plot
        shrink = 0.75
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * shrink, box.height])
        
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * shrink, box.height])
        
        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * shrink, box.height])
        
        box = ax4.get_position()
        ax4.set_position([box.x0, box.y0, box.width * shrink, box.height])
        
        #___Energy Plots ___________________________________________________________________
        wts = np.linspace(4,2,3)
        ax1.plot(t_array,s_array.m,    linewidth = wts[0], label='TES Reserve Quantity')
        ax1.plot(t_array,ucsu_array.m, linewidth = wts[1], label='Cycle Startup Energy Inventory')
        ax1.plot(t_array,ursu_array.m, linewidth = wts[2], label='Nuclear Startup Energy Inventory')
        
        #-Line marking the 24 hour line
        ax1.plot(time_24hr, energy_24hr, 'k--', linewidth=lw )
        
        #-Labels
        ax1.set_ylabel('Energy (MWh)', labelpad=lp, fontsize=fs, fontweight='bold')
        ax1.legend(loc=loc,fontsize=fsl,bbox_to_anchor=(1.4, 1.0)) #bbox stuff places legend off-plot
        
        #___Power Plots ___________________________________________________________________
        wts = np.linspace(6,2,6)
        ax2.plot(t_array,wdot_array.m, linewidth = wts[0], label='Cycle Out (E)')
        ax2.plot(t_array,x_array.m,    linewidth = wts[1], label='Cycle In (T)')
        ax2.plot(t_array,xr_array.m,   linewidth = wts[2], label='Nuclear Out (T)')
        ax2.plot(t_array,xrsu_array.m, linewidth = wts[3], label='Nuclear Startup (T)')
        ax2.plot(t_array,wdot_s_array.m, linewidth = wts[4], label='Energy Sold to Grid (E)')
        ax2.plot(t_array,wdot_p_array.m, linewidth = wts[5], label='Energy Purchased (E)')
        
        #-Line marking the 24 hour line
        ax2.plot(time_24hr, power_24hr, 'k--', linewidth=lw )
        
        #-Labels
        ax2.set_ylabel('Power (MW)', labelpad=lp, fontsize=fs, fontweight='bold')
        ax2.legend(loc=loc,fontsize=fsl,bbox_to_anchor=(1.42, 1.0))
        
        #___Power Plots x2 ___________________________________________________________________
        wts = np.linspace(6,2,6)
        ax3.plot(t_array,wdot_delta_plus_array.m,  linewidth = wts[0], label='PC Ramp Up (E)')
        ax3.plot(t_array,wdot_delta_minus_array.m, linewidth = wts[1], label='PC Ramp Down (E)')
        ax3.plot(t_array,wdot_v_plus_array.m,      linewidth = wts[2], label='PC Ramp Up Beyond (E')
        ax3.plot(t_array,wdot_v_minus_array.m,     linewidth = wts[3], label='PC Ramp Down Beyond(E')
        # ax3.plot(t_array,wdot_s_prev_delta_plus_array.m,  linewidth = wts[4], label='UB on Delta w (E')
        # ax3.plot(t_array,wdot_s_prev_delta_minus_array.m, linewidth = wts[5], label='LB on Delta w (E')
        
        #-Line marking the 24 hour line
        ax3.plot(time_24hr, power_24hr, 'k--', linewidth=lw )
        
        #-Labels
        ax3.set_ylabel('Power (MW)', labelpad=lp, fontsize=fs, fontweight='bold')
        ax3.legend(loc=loc,fontsize=fsl,bbox_to_anchor=(1.34, 1.0))
        
        #___Price Plot ________________________________________________________________
        ax5.bar(t_array, p_array, np.diff(t_array)[0], alpha=0.35, label="Price")
        
        wts = np.linspace(10,1,14)
        #___Binary Plots ______________________________________________________________
        ax4.plot(t_array, yr_array,    linewidth = wts[0], label='Receiver On?')
        ax4.plot(t_array, yrhsp_array, linewidth = wts[1], label='Receiver HSU Pen?')
        ax4.plot(t_array, yrsb_array,  linewidth = wts[2], label='Receiver SB?')
        ax4.plot(t_array, yrsd_array,  linewidth = wts[3], label='Receiver SD?')
        ax4.plot(t_array, yrsu_array,  linewidth = wts[4], label='Receiver SU?')
        ax4.plot(t_array, yrsup_array, linewidth = wts[5], label='Receiver CSU Pen?')
        ax4.plot(t_array, y_array,     linewidth = wts[6], label='Cycle Generating Power?')
        ax4.plot(t_array, ychsp_array, linewidth = wts[7], label='Cycle HSU Pen?')
        ax4.plot(t_array, ycsb_array,  linewidth = wts[8], label='Cycle SB?')
        ax4.plot(t_array, ycsd_array,  linewidth = wts[9], label='Cycle SD?')
        ax4.plot(t_array, ycsu_array,  color='k', linewidth = wts[10], label='Cycle SU?')
        ax4.plot(t_array, ycsup_array, linewidth = wts[11], label='Cycle CSU Pen?')
        ax4.plot(t_array, ycgb_array,  linewidth = wts[12], label='Cycle Began Gen?')
        ax4.plot(t_array, ycge_array,  linewidth = wts[13], label='Cycle Stopped Gen?')
        ax4.plot(t_array, binary_div, 'k--', linewidth=2)
        
        #-Line marking the 24 hour line
        ax4.plot(time_24hr, binary_24hr, 'k--', linewidth=lw )
        
        #-Labels
        ax4.set_xlabel('Time (hrs)', labelpad=lp, fontsize=fs, fontweight='bold')
        ax4.set_ylabel('Binary Vars', labelpad=lp, fontsize=fs, fontweight='bold')
        ax4.legend(loc=loc,fontsize=fsl,bbox_to_anchor=(1.35, 1.3) )
        ax5.set_ylabel('Price ($/kWh)', labelpad=lps, fontsize=fs, fontweight='bold')
        ax5.legend(loc=loc_cr,fontsize=fsl)
        
        #-Axis limits and tickmarks
        ax4.set_ylim(-0.2,3.4)
        ax4.set_yticks([0,1,2,3])
        ax4.set_yticklabels([0,1,0,1])
        ax5.set_ylim(0.4,2.5)
        ax5.set_yticks([0.50,0.75, 1])
        
        
        # set full title 
        ax1.set_title('Pyomo Results - Normal',fontweight='bold')
        
        
        # plt.tight_layout()
        
        
        # =============================================================================
        # objective function
        # =============================================================================
        
        # fig = plt.figure()
        # ax1o = fig.add_subplot(211)
        # # ax1o = fig.add_subplot(111)
        # ax2o = fig.add_subplot(212)
        
        # #profit
        # Obj_profit      = 0.1 * d_array * p_array * ( wdot_s_array - wdot_p_array ) 
        # Obj_cycle_susd  = (Ccsu * ycsup_array + 0.1*Cchsp*ychsp_array + ycsd_array)
        # Obj_cycle_ramp  = C_delta_w*(wdot_delta_plus_array + wdot_delta_minus_array) + C_v_w*(wdot_v_plus_array + wdot_v_minus_array)
        # Obj_rec_susd    = Crsu*yrsup_array + Crhsp*yrhsp_array + yrsd_array
        # Obj_ops         = d_array*(Cpc*wdot_array + Ccsb*Qb*ycsb_array + Crec*xr_array)
        
        # ax1o.plot(t_array, Obj_profit, linewidth=lw, label='Profit Term')
        # ax2o.plot(t_array, -Obj_cycle_susd, linewidth=lw, label='Cycle Startup/Shutdown Term')
        # ax1o.plot(t_array, -Obj_cycle_ramp, linewidth=lw, label='Cycle Ramping Term')
        # ax1o.plot(t_array, -Obj_rec_susd, linewidth=lw, label='Receiver Startup/Shutdown Term')
        # ax1o.plot(t_array, -Obj_ops, linewidth=lw, label='Cycle and Rec Ops Term')
        
        # ax1o.legend(loc='best')
        # ax2o.legend(loc='best')
        
        # ax1o.set_title('Pyomo Obj Fun - Normal', fontweight='bold')
        
        # plt.tight_layout()

        
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

    