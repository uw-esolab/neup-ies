#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 14:46:10 2022

@author: gabrielsoto
"""

from pylab import rc
import matplotlib.pyplot as plt
import pyomo.environ as pe
import numpy as np
import pint, copy
u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
rc('axes', linewidth=2)
rc('font', weight='bold', size=12)
from util.PostProcessing import OutputExtraction
from util.PostProcessing import Plots
from util.PostProcessing import DispatchPlots
from util.SolarPlots import SolarOutputExtraction
from util.SolarPlots import SolarPlots
from util.SolarPlots import SolarDispatchPlots

class DualOutputExtraction(SolarOutputExtraction):
    """
    The OutputExtraction class is a part of the PostProcessing family of classes.
    It extracts outputs from SSC or Dispatch models. Can be called from other 
    plotting classes to extract outputs.  
    """
    
    def __init__(self, module):
        """ Initializes the OutputExtraction module

        The instantiation of this class receives a full module object, the module
        being one of the NE2 modules in the /neup-ies/simulations/modules directory.
        or a Dispatch module. 

        Inputs:
            module (object)     : object representing NE2 module class after simulations
        """
        
        SolarOutputExtraction.__init__(self, module)


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
        self.q_thermal    = np.asarray(mod_out.Q_thermal) * u.MW
        self.q_dot_nuc_in = np.asarray(mod_out.q_dot_nuc_inc) * u.MW
        self.q_nuc_thermal    = np.asarray(mod_out.Q_nuc_thermal) * u.MW
        self.q_pb         = np.asarray(mod_out.q_pb) * u.MW
        self.q_dot_pc_su  = np.asarray(mod_out.q_dot_pc_startup) * u.MW
        self.m_dot_pc     = np.asarray(mod_out.m_dot_pc) * u.kg/u.s
        self.m_dot_rec    = np.asarray(mod_out.m_dot_rec) * u.kg/u.s
        self.m_dot_nuc    = np.asarray(mod_out.m_dot_nuc) * u.kg/u.s
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
        self.q_nuc_design  = self.mod.SystemDesign.q_dot_nuclear_des * u.MW  # receiver design thermal power
        self.p_pb_design   = self.mod.SystemDesign.P_ref * u.MW              # power block design electrical power
        self.eta_design    = self.mod.SystemDesign.design_eff                # power block design efficiency
        self.q_pb_design  = (self.p_pb_design / self.eta_design).to('MW')    # power block design thermal rating
        self.T_htf_hot    = (self.mod.SystemDesign.T_htf_hot_des*u.celsius).to('degK')   # heat transfer fluid Hot temp
        self.T_htf_cold   = (self.mod.SystemDesign.T_htf_cold_des*u.celsius).to('degK')  # heat transfer fluid Cold temp
        self.e_tes_design = (self.q_pb_design * self.mod.SystemDesign.tshours*u.hr).to('MWh')  # TES storage capacity (kWht)

        # operating modes
        op_mode_result, modes_order = np.unique(self.op_mode_1, return_index=True) # mode orders and re-ordering
        self.op_mode_result = op_mode_result[np.argsort(modes_order)]     # re-order modes by first appearance of each


    def set_pyomo_outputs(self):
        """ Method to define list of Pyomo output arrays

        This method extracts outputs from the Pyomo Dispatch model, converts them to numpy arrays
        and saves them to `self`. 
        """
        
        OutputExtraction.set_pyomo_outputs(self)
        SolarOutputExtraction.set_pyomo_outputs(self)


class DualPlots(SolarPlots):
    """
    The Plots class is a part of the PostProcessing family of classes. It can 
    output results from SSCmodels. This might be a 
    Sisyphus-ian task as some of the parameters are very unique to whatever 
    plot you are trying to make, but will do my best to provide basic structures
    and templates to work off of. 

    Note that the Plots class must be initialized before using. 
    """

    def __init__(self, module, **kwargs):
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
    
        SolarPlots.__init__(self, module, **kwargs)
        
        
    def set_extractor(self):
        """ Setting the output extraction class
        """
        
        self.extractor = DualOutputExtraction
        

    def plot_SSC_power_and_energy(self, ax=None, title_label=None, plot_all_time=True, \
                                  start_hr=0, end_hr=48, hide_x=False, x_legend=1.2, \
                                  y_legend_L=1.0, y_legend_R=1.0):
        """ Method to plot power and energy data on single plot

        This method is used specifically to plot power and energy data from SSC simulation
        results. Built-in options to plot legend off-axis. 

        Inputs:
            ax (object)         : axis object to plot on
            plot_all_time(bool) : are we plotting all results or just a portion?
            title_label(str)    : title name for plot
            start_hr (int)      : (plot_all_time==False) hour used for starting index 
            end_hr (int)        : (plot_all_time==False) hour used for ending index
            hide_x(bool)        : hiding the x-axis from this particular plot
            x_legend (float)    : (legend_offset==True) x-offset defining left-most side of legend
            y_legend_L (float)  : (legend_offset==True) y-offset of left y-axis plot
            y_legend_R (float)  : (legend_offset==True) y-offset of right y-axis plot
        """
        #========================#
        #--- Creating Figure  ---#
        #========================#

        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10, 5])
            ax = fig.gca()   # this is the power plot

        # twin axis to plot energy on opposite y-axis
        ax2 = ax.twinx()  # this is the energy plot
        
        # custom y limits and ticks to be integers for Power
        #TODO: add these as some sort of input, maybe dict?
        ax.set_ylim(-100, 1100)
        ax.set_yticks([0, 250, 500, 750, 1000])

        # plot Power arrays
        power_array_list = ['p_cycle', 'q_dot_rec_in', 'q_thermal', 'q_nuc_thermal', 'gen', 'q_dot_pc_su'] # list of array strings
        power_label_list = ['P_cycle (Electric)',
                            'Q_dot Receiver Incident (Thermal)',
                            'Q_dot to Salt from CSP (Thermal)',
                            'Q_dot to Salt from LFR (Thermal)',
                            'Power generated (Electric)',
                            'PC startup thermal power (Thermal)'] # list of labels for each array string to extract from Outputs
        power_ylabel = 'Power \n(MW)'
        ax  = self.plot_SSC_generic(ax, array_list=power_array_list, \
                                    label_list=power_label_list, \
                                    y_label=power_ylabel, \
                                    title_label=title_label, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x)

        # custom y limits and ticks to be integers for Energy
        ax2.set_ylim(-0.05*self.e_tes_design.m, 1.05*self.e_tes_design.m)

        # plot Energy array(s)
        energy_array_list = ['e_ch_tes']
        energy_label_list = ['Salt Charge Energy Level (Thermal)']
        energy_ylabel = 'Energy \n(MWh)'
        ax2 = self.plot_SSC_generic(ax2, array_list=energy_array_list, \
                                    label_list=energy_label_list, \
                                    y_label=energy_ylabel, \
                                    title_label=None, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x, left_axis=False)
        
        # set line color to default C4 (purple)
        ax2.get_lines()[0].set_color("C6")
        
        #========================#
        #---- Setting Labels ----#
        #========================#
        
        # customizing legend(s)
        if self.legend_offset:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl, bbox_to_anchor=(x_legend, y_legend_L)) # plot legend for Power arrays
            ax2.legend(loc=self.loc_ul, fontsize=self.fsl, bbox_to_anchor=(x_legend, y_legend_R)) # plot legend for Energy arrays and also 
        else:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl) # plot legend for Power arrays
            ax2.legend(loc=self.loc_ul, fontsize=self.fsl) # plot legend for Energy arrays and also 


    def plot_SSC_op_modes(self, ax=None, title_label=None, plot_all_time=True, \
                          start_hr=0, end_hr=48, hide_x=False,  x_legend=1.2, \
                          y_legend_L=1.0, y_legend_R=1.0):
        """ Method to plot operating modes history on single plot

        This method is used specifically to plot operating modes and relative 
        pricing data from SSC simulation results. Built-in options to plot legend off-axis. 

        Inputs:
            ax (object)         : axis object to plot on
            plot_all_time(bool) : are we plotting all results or just a portion?
            title_label(str)    : title name for plot
            start_hr (int)      : (plot_all_time==False) hour used for starting index 
            end_hr (int)        : (plot_all_time==False) hour used for ending index
            hide_x(bool)        : hiding the x-axis from this particular plot
            x_legend (float)    : (legend_offset==True) x-offset defining left-most side of legend
            y_legend_L (float)  : (legend_offset==True) y-offset of left y-axis plot
            y_legend_R (float)  : (legend_offset==True) y-offset of right y-axis plot
        """
        #========================#
        #--- Creating Figure  ---#
        #========================#

        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10, 5])
            ax = fig.gca()   # this is the Op Modes plot

        # twin axis to plot pricing on opposite y-axis
        axO = ax  # this is the OP modes plot
        ax2 = axO.twinx()  # this is the pricing plot

        # custom y limits and ticks to be integers
        price = self.price[ start_hr:end_hr ]
        minP = 0
        maxP = 1.05*price.max()
        spacing = 0.5 if maxP-minP < 2 else 1
        
        ax2.set_ylim(minP, maxP)
        ax2.set_yticks( np.arange( minP, maxP, spacing) )
        
        # plot price array(s)
        price_array_list = ['price']
        price_label_list = [None]
        price_ylabel = 'Tariff \n($/kWh)'
        ax2 = self.plot_SSC_generic(ax2, array_list=price_array_list, \
                                    label_list=price_label_list, \
                                    y_label=price_ylabel, \
                                    title_label=None, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x, \
                                    is_bar_graph=True, left_axis=False)
        
        #========================#
        #--- Rearranging Modes --#
        #========================#
        
        # time series of operating mode indeces
        op_mode_1 = self.op_mode_1
        
        # array of unique mode indeces and the time-order in which they appear
        op_mode_unique, op_mode_order = np.unique(op_mode_1, return_index=True) 
        
        # arranging unique modes by order in which they appear
        op_mode_unique = op_mode_unique[np.argsort(op_mode_order)]   
        
        new_order = np.arange(1, len(op_mode_unique)+1)
        op_mode_ordered = copy.deepcopy( op_mode_1 )
        for i,op in enumerate(op_mode_unique):
            op_mode_ordered[np.where(op_mode_ordered == op)] = new_order[i]
        
        
        self.op_mode_ordered = op_mode_ordered

        
        #========================#
        #--- original plot     --#
        #========================#

        # plot operating mode arrays
        op_array_list = ['op_mode_ordered'] # list of array strings
        op_label_list = [None]        # list of labels for each array string to extract from Outputs
        op_ylabel = '' # 'Operating \nMode'
        axO, d_slice, t_plot  = self.plot_SSC_generic(axO, array_list=op_array_list, \
                                    label_list=op_label_list, \
                                    y_label=op_ylabel, \
                                    title_label=title_label, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x, \
                                    return_extra=True)
            
        # custom y limits and ticks to be integers
        axO.set_ylim(0, len(op_mode_unique)+2 )
        axO.set_yticks(np.arange(0, len(op_mode_unique)+2, 1))
        
        #==================================================#
        #---- extract operating modes to designated array
        #==================================================#
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        colors = np.tile( colors, 12 )
        
        # time series of operating mode indeces
        op_mode_1 = self.op_mode_1[d_slice]
        
        # array of unique mode indeces and the time-order in which they appear
        op_mode_unique, op_mode_order = np.unique(op_mode_1, return_index=True) 
        
        # arranging unique modes by order in which they appear
        op_mode_unique = op_mode_unique[np.argsort(op_mode_order)]   
        
        # Plotting data points over the OP mode line with different colors and labels
        count = 0
        for op in op_mode_unique:
            # getting unique operating modes
            inds = (op_mode_1 == op)
            # individual index getting plotted with unique color and label
            if np.sum(inds):
                axO.plot(t_plot[inds].m, (count+1)*np.ones(np.sum(inds)), 'o', color=colors[count], label=self.operating_modes[int(op)])
                count +=1
                
        #==================================================#
        #---- handling labels
        #==================================================#
        
        # list of unique op mode strings
        op_modes_labels = [self.operating_modes[int(l)] for l in op_mode_unique]
        
        # split strings by "__" and create dictionary of all subsystems (PC, TES, etc)
        op_modes_subsystem_splits = [op.split('__') for op in op_modes_labels]
        op_modes_dict = {subsys.split('_')[0]:[] for subsys in op_modes_subsystem_splits[0]}
        
        for op_mode in op_modes_subsystem_splits:
            for subsys in op_mode:
                op_modes_dict[subsys.split('_')[0]].append( subsys.split('_')[1] )
        
        max_chars = np.array([ len( max(op_modes_dict[subsys], key=len) ) for subsys in op_modes_dict.keys()])
        new_tick_labels = []
        for label in op_modes_labels:
            formatted_label = ''
            for i,(c,s) in enumerate( zip(max_chars, op_modes_dict.keys()) ):
                subsys_mode = label.split('__')[i].split('_')[1:]
                subsys_mode_label = '_'.join(subsys_mode)
                subsys_mode_label = subsys_mode_label.ljust(c+2,'_')
                
                # subsys_op_mode_label = "{0}_{1}".format(s,subsys_mode_label )
                formatted_label += "{0}_{1}".format(s,subsys_mode_label ) + '__'
            new_tick_labels.append( formatted_label )
        
        blahx = axO.set_yticks( np.arange(1,count+1) )
        axO.set_yticklabels(new_tick_labels)
        for i,(ytick,color) in enumerate( zip(axO.get_yticklabels(), colors) ):
            
            current_label = "OFF"
            
            count = 0
            for j, _ in enumerate(new_tick_labels[i]):
                if new_tick_labels[i][j:j + len('OFF')] == 'OFF':
                    count +=1

            if count == 3:
                ytick.set_color('k')
            else:
                ytick.set_color(color)

        #========================#
        #---- Setting Labels ----#
        #========================#
        
        # set OP Mode line color to black
        axO.get_lines()[0].set_color("w")
        
        plt.tight_layout()
        
        # customizing legend(s)
        # if self.legend_offset:
        #     ax.legend( loc=self.loc_ul, fontsize=self.fsl, bbox_to_anchor=(x_legend, y_legend_L)) # plot legend for Op Modes arrays
        # else:
        #     ax.legend( loc=self.loc_ul, fontsize=self.fsl) # plot legend for Op Modes arrays



    def set_operating_modes_list(self):
        """ Method to define list of operating modes

        This method creates a list of operating modes pertaining to the specific SSC module we are
        using. This particular list was taken from /ssc/tcs/csp_solver_core.h, in the 
        `C_system_operating_modes` class.
        """

        self.operating_modes = [
            'ITER_START',
            'CR_OFF__PC_OFF__TES_OFF__NUC_ON',
            'CR_SU__PC_OFF__TES_OFF__NUC_ON',
            'CR_ON__PC_SU__TES_OFF__NUC_ON',
            'CR_ON__PC_SB__TES_OFF__NUC_ON',
            'CR_ON__PC_RM_HI__TES_OFF__NUC_ON',
            'CR_ON__PC_RM_LO__TES_OFF__NUC_ON',
            'CR_DF__PC_MAX__TES_OFF__NUC_ON',
            'CR_OFF__PC_SU__TES_DC__NUC_ON',
            'CR_ON__PC_OFF__TES_CH__NUC_ON',
            'CR_ON__PC_TARGET__TES_CH__NUC_ON',
            'CR_ON__PC_TARGET__TES_DC__NUC_ON',
            'CR_ON__PC_RM_LO__TES_EMPTY__NUC_ON',
            'CR_DF__PC_OFF__TES_FULL__NUC_ON',
            'CR_OFF__PC_SB__TES_DC__NUC_ON',
            'CR_OFF__PC_MIN__TES_EMPTY__NUC_ON',
            'CR_OFF__PC_RM_LO__TES_EMPTY__NUC_ON',
            'CR_ON__PC_SB__TES_CH__NUC_ON',
            'CR_SU__PC_MIN__TES_EMPTY__NUC_ON',
            'CR_SU__PC_SB__TES_DC__NUC_ON',
            'CR_ON__PC_SB__TES_DC__NUC_ON',
            'CR_OFF__PC_TARGET__TES_DC__NUC_ON',
            'CR_SU__PC_TARGET__TES_DC__NUC_ON',
            'CR_ON__PC_RM_HI__TES_FULL__NUC_ON',
            'CR_ON__PC_MIN__TES_EMPTY__NUC_ON',
            'CR_SU__PC_RM_LO__TES_EMPTY__NUC_ON',
            'CR_DF__PC_MAX__TES_FULL__NUC_ON',
            'CR_ON__PC_SB__TES_FULL__NUC_ON',
            'CR_SU__PC_SU__TES_DC__NUC_ON',
            'CR_ON__PC_SU__TES_CH__NUC_ON',
            'CR_DF__PC_SU__TES_FULL__NUC_ON',
            'CR_DF__PC_SU__TES_OFF__NUC_ON',
            'CR_TO_COLD__PC_TARGET__TES_DC__NUC_ON',
            'CR_TO_COLD__PC_RM_LO__TES_EMPTY__NUC_ON',
            'CR_TO_COLD__PC_SB__TES_DC__NUC_ON',
            'CR_TO_COLD__PC_MIN__TES_EMPTY__NUC_ON',
            'CR_TO_COLD__PC_OFF__TES_OFF__NUC_ON',
            'CR_TO_COLD__PC_SU__TES_DC__NUC_ON',
            'CR_OFF__PC_OFF__TES_CH__NUC_ON',
            'CR_OFF__PC_SU__TES_CH__NUC_ON',
            'CR_OFF__PC_SU__TES_OFF__NUC_ON',
            'CR_SU__PC_OFF__TES_CH__NUC_ON',
            'CR_SU__PC_SU__TES_CH__NUC_ON',
            'CR_SU__PC_SU__TES_OFF__NUC_ON',
            'CR_OFF__PC_SB__TES_OFF__NUC_ON',
            'CR_OFF__PC_SB__TES_CH__NUC_ON',
            'CR_OFF__PC_SB__TES_FULL__NUC_ON',
            'CR_OFF__PC_RM_LO__TES_OFF__NUC_ON',
            'CR_OFF__PC_TARGET__TES_CH__NUC_ON',
            'CR_OFF__PC_RM_HI__TES_OFF__NUC_ON',
            'CR_OFF__PC_RM_HI__TES_FULL__NUC_ON',
            'CR_SU__PC_SB__TES_OFF__NUC_ON',
            'CR_SU__PC_SB__TES_CH__NUC_ON',
            'CR_SU__PC_SB__TES_FULL__NUC_ON',
            'CR_SU__PC_RM_LO__TES_OFF__NUC_ON',
            'CR_SU__PC_TARGET__TES_CH__NUC_ON',
            'CR_SU__PC_RM_HI__TES_OFF__NUC_ON',
            'CR_SU__PC_RM_HI__TES_FULL__NUC_ON',
            'ITER_END']


class DualDispatchPlots(DualPlots):
    """
    The Plots class is a part of the PostProcessing family of classes. It can 
    output results from Pyomo Dispatch models. 

    Note that the DispatchPlots class must be initialized before using. 
    """

    def __init__(self, module, **kwargs):
        """ Initializes the Plots module

        The instantiation of this class receives a full Dispatch object, the module
        being one of the Pyomo Dispatch model created in the /neup-ies/simulations/dispatch directory.
        It also contains various inputs relating to cosmetic parameters for
        matplotlib plots. 

        Inputs:
            module (object)     : object representing Pyomo Dispatch Model with results
            fsl (str)           : fontsize for legend
            loc (str)           : location of legend
            legend_offset(bool) : are we plotting legends off-axis?
            lp (int)            : labelpad for axis labels
            lps (int)           : labelpad for axis labels - short version
            fs (int)            : fontsize for labels, titles, etc.
            lw (int)            : linewidth for plotting
            x_shrink (float)    : (legend_offset==True) amount to shrink axis to make room for legend
        """
        
        # initialize Plots class
        SolarDispatchPlots.__init__( self, module, **kwargs)
        
    
    def set_plotter(self):
        """ Setting class for plotting
        """
        
        self.plotter = DualPlots


    def plot_pyomo_energy(self, ax=None, **kwargs):
        """ Method to plot energy data from Pyomo Dispatch on single plot

        This method is used specifically to plot energy data from Dispatch simulation
        results. Built-in options to plot legend off-axis. 

        Inputs:
            ax (object)         : axis object to plot on
            plot_all_time(bool) : are we plotting all results or just a portion?
            title_label(str)    : title name for plot
            start_hr (int)      : (plot_all_time==False) hour used for starting index 
            end_hr (int)        : (plot_all_time==False) hour used for ending index
            hide_x(bool)        : hiding the x-axis from this particular plot
            x_legend (float)    : (legend_offset==True) x-offset defining left-most side of legend
            y_legend_L (float)  : (legend_offset==True) y-offset of left y-axis plot
            y_legend_R (float)  : (legend_offset==True) y-offset of right y-axis plot
        """
        
        SolarDispatchPlots.plot_pyomo_energy(self, ax=ax, **kwargs)
        

    def plot_pyomo_power(self, ax=None, title_label=None, plot_all_time=True, \
                          start_hr=0, end_hr=48, hide_x=False,  x_legend=1.2, \
                          y_legend_L=1.0, y_legend_R=1.0):
        """ Method to plot power and pricing data from Pyomo Dispatch on single plot

        This method is used specifically to plot power data from Dispatch simulation
        results. Built-in options to plot legend off-axis. 

        Inputs:
            ax (object)         : axis object to plot on
            plot_all_time(bool) : are we plotting all results or just a portion?
            title_label(str)    : title name for plot
            start_hr (int)      : (plot_all_time==False) hour used for starting index 
            end_hr (int)        : (plot_all_time==False) hour used for ending index
            hide_x(bool)        : hiding the x-axis from this particular plot
            x_legend (float)    : (legend_offset==True) x-offset defining left-most side of legend
            y_legend_L (float)  : (legend_offset==True) y-offset of left y-axis plot
            y_legend_R (float)  : (legend_offset==True) y-offset of right y-axis plot
        """

        #========================#
        #--- Creating Figure  ---#
        #========================#

        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10, 5])
            ax = fig.gca()   # this is the power plot
    
        # twin axis to plot pricing on opposite y-axis
        ax2 = ax.twinx()  # this is the pricing plot

        # plot energy arrays
        power_array_list = ['wdot_array', 'x_array', 'xr_array',
                            'xrsu_array', 'xn_array',
                            'wdot_s_array', 'wdot_p_array'] # list of array strings
        power_label_list = ['Cycle Out (E)', 'Cycle In (T)',
                            'Receiver Out (T)', 'Receiver Startup (T)', 'Nuclear Out (T)',
                            'Energy Sold to Grid (E)', 'Energy Purchased (E)'] # list of labels for each array string to extract from Outputs
        power_wts = np.linspace(6, 2, 6).tolist()
        power_ylabel = 'Power \n(MW)'

        ax  = self.plot_SSC_generic(ax, array_list=power_array_list, \
                                    label_list=power_label_list, \
                                    y_label=power_ylabel, \
                                    lw_list=power_wts, \
                                    title_label=title_label, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x)
        
        #========================#
        #--- Extra arrays  ------#
        #========================#
        
        # retrieving arrays
        wdot_array    = self.wdot_array.m    
        x_array       = self.x_array.m
        xr_array      = self.xr_array.m
        xrsu_array    = self.xrsu_array.m    
        wdot_s_array  = self.wdot_s_array.m
        wdot_p_array  = self.wdot_p_array.m
        
        # vertical power line at midpoint
        power_vert = np.linspace(  np.min([wdot_array, x_array, xr_array, \
                                           xrsu_array, wdot_s_array, wdot_p_array]), 
                                   np.max([wdot_array, x_array, xr_array, \
                                           xrsu_array, wdot_s_array, wdot_p_array])*1.1, 
                                   self.T )
            
        # Line marking the midpoint line
        ax.plot(self.time_midway, power_vert, 'k--', linewidth=self.lw)   
            
        #========================#
        #--- Double axis      ---#
        #========================#
        
        # plot price array(s)
        price_array_list = ['p_array']
        price_label_list = [None]
        price_ylabel = 'Tariff \n($/kWh)'
        ax2 = self.plot_SSC_generic(ax2, array_list=price_array_list, \
                                    label_list=price_label_list, \
                                    y_label=price_ylabel, \
                                    title_label=None, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x, \
                                    is_bar_graph=True, left_axis=False)
            
        #========================#
        #---- Setting Labels ----#
        #========================#
        
        # customizing legend(s)
        if self.legend_offset:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl, bbox_to_anchor=(x_legend, y_legend_L)) # plot legend for Power arrays
        else:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl) # plot legend for Power arrays
            

    def plot_pyomo_dual_bin(self, ax=None, title_label=None, plot_all_time=True, \
                          start_hr=0, end_hr=48, hide_x=False,  x_legend=1.2, \
                          y_legend_L=1.0, y_legend_R=1.0):
        """ Method to plot solar binary data from Pyomo Dispatch on single plot

        This method is used specifically to plot data from Dispatch simulation
        results pertaining to nuclear binary values. Variables include whether
        nuclear plant is running (0 or 1), etc. Built-in options to plot legend 
        off-axis. 

        Inputs:
            ax (object)         : axis object to plot on
            plot_all_time(bool) : are we plotting all results or just a portion?
            title_label(str)    : title name for plot
            start_hr (int)      : (plot_all_time==False) hour used for starting index 
            end_hr (int)        : (plot_all_time==False) hour used for ending index
            hide_x(bool)        : hiding the x-axis from this particular plot
            x_legend (float)    : (legend_offset==True) x-offset defining left-most side of legend
            y_legend_L (float)  : (legend_offset==True) y-offset of left y-axis plot
            y_legend_R (float)  : (legend_offset==True) y-offset of right y-axis plot
        """

        #========================#
        #--- Creating Figure  ---#
        #========================#

        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10, 5])
            ax = fig.gca()   # this is the power plot

        # plot nuclear binary arrays
        recbin_array_list = ['yn_array', 'yr_array',   'yrhsp_array', 'yrsb_array',
                             'yrsd_array', 'yrsu_array',  'yrsup_array' ] # list of array strings
        recbin_label_list = ['Is Nuclear On?', 'Is Receiver On?', 'Is Receiver HSU Pen?', 
                             'Is Receiver SB?', 'Is Receiver SD?', 
                             'Is Receiver SU?', 'Is Receiver CSU Pen?' ] # list of labels for each array string to extract from Outputs
        recbin_wts = np.linspace(10, 1.5, 6).tolist()
        recbin_ylabel = 'CSP and LFR \nBinary \nVariables'

        ax  = self.plot_SSC_generic(ax, array_list=recbin_array_list, \
                                    label_list=recbin_label_list, \
                                    y_label=recbin_ylabel, \
                                    lw_list=recbin_wts, \
                                    title_label=title_label, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x)
            
        #========================#
        #--- Extra arrays  ---#
        #========================#

        # vertical binary line at midpoint
        binary_vert = np.linspace( -0.5, 1.5, self.T )
            
        # Line marking the midpoint line
        ax.plot(self.time_midway, binary_vert, 'k--', linewidth=self.lw)
        
        # set binary ticks and labels
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['No', 'Yes'])

        #========================#
        #---- Setting Labels ----#
        #========================#
        
        # customizing legend(s)
        if self.legend_offset:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl, bbox_to_anchor=(x_legend, y_legend_L)) # plot legend for Nuclear arrays
        else:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl) # plot legend for Nuclear arrays


    def plot_pyomo_cycle_bin(self, ax=None, **kwargs):
        """ Method to plot power cycle binary data from Pyomo Dispatch on single plot

        This method is used specifically to plot data from Dispatch simulation
        results pertaining to cycle binary values. Variables include whether
        power cycle is running (0 or 1), etc. Built-in options to plot legend 
        off-axis. 

        Inputs:
            ax (object)         : axis object to plot on
            plot_all_time(bool) : are we plotting all results or just a portion?
            title_label(str)    : title name for plot
            start_hr (int)      : (plot_all_time==False) hour used for starting index 
            end_hr (int)        : (plot_all_time==False) hour used for ending index
            hide_x(bool)        : hiding the x-axis from this particular plot
            x_legend (float)    : (legend_offset==True) x-offset defining left-most side of legend
            y_legend_L (float)  : (legend_offset==True) y-offset of left y-axis plot
            y_legend_R (float)  : (legend_offset==True) y-offset of right y-axis plot
        """
        
        DispatchPlots.plot_pyomo_cycle_bin(self, ax=ax, **kwargs)