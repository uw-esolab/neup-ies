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
from util.PostProcessing import OutputExtraction
from util.PostProcessing import Plots
from util.PostProcessing import DispatchPlots


class SolarOutputExtraction(OutputExtraction):
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
        
        OutputExtraction.__init__(self, module)
        
        
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


    def set_pyomo_outputs(self):
        """ Method to define list of Pyomo output arrays

        This method extracts outputs from the Pyomo Dispatch model, converts them to numpy arrays
        and saves them to `self`. 
        """
        
        # Dispatch Model with results
        dm = self.dm
        u = self.u
        
        # lambda functions
        extract_from_model = lambda name: getattr(dm.model, name)
        extract_array      = lambda name: np.array([ pe.value(extract_from_model(name)[t]) for t in dm.model.T ])
        extract_energy     = lambda name: (extract_array(name)*u.kWh).to('MWh') 
        extract_power      = lambda name: (extract_array(name)*u.kW).to('MW') 
        
        # Time and Pricing Arrays
        self.t_full = extract_array('Delta_e') * u.hr
        self.p_array = extract_array('P')
        
        # marking the midway point
        self.T = len(self.t_full)
        self.time_midway = np.ones([self.T])*int(self.T/2)

        # Energy Arrays
        self.s_array     = extract_energy('s')
        self.ucsu_array  = extract_energy('ucsu')
        self.ursu_array  = extract_energy('ursu')
        
        # Power Arrays
        self.wdot_array             = extract_power('wdot')
        self.wdot_delta_plus_array  = extract_power('wdot_delta_plus')
        self.wdot_delta_minus_array = extract_power('wdot_delta_minus')
        self.wdot_v_plus_array      = extract_power('wdot_v_plus')
        self.wdot_v_minus_array     = extract_power('wdot_v_minus')
        self.wdot_s_array           = extract_power('wdot_s')
        self.wdot_p_array           = extract_power('wdot_p')
        self.x_array                = extract_power('x')
        self.xr_array               = extract_power('xr')
        self.xrsu_array             = extract_power('xrsu')

        # Binary Arrays (Nuclear)
        self.yr_array    = extract_array('yr') 
        self.yrhsp_array = extract_array('yrhsp') 
        self.yrsb_array  = extract_array('yrsb') 
        self.yrsd_array  = extract_array('yrsu') 
        self.yrsu_array  = extract_array('yrsu') 
        self.yrsup_array = extract_array('yrsup') 
        
        # Binary Arrays (Cycle)
        self.y_array     = extract_array('y') 
        self.ychsp_array = extract_array('ychsp') 
        self.ycsb_array  = extract_array('ycsb') 
        self.ycsd_array  = extract_array('ycsd') 
        self.ycsu_array  = extract_array('ycsu') 
        self.ycsup_array = extract_array('ycsup') 
        self.ycgb_array  = extract_array('ycgb')  
        self.ycge_array  = extract_array('ycge') 


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
        
        self.u = u

        # user-defined plotting parameters
        self.lp  = lp   # labelpad
        self.lps = lps  # labelpad short
        self.fs  = fs   # fontsize
        self.lw  = lw   # linewidth
        self.fsl = fsl  # fontsize legend
        self.loc = loc  # location of legend

        # offsetting legend
        self.legend_offset = legend_offset # boolean - are we plotting legends off-axis?
        self.x_shrink      = x_shrink # amount to shrink x-asis by to make room for legend

        # alternate legend locations
        self.loc_ur = 'upper right'
        self.loc_ul = 'upper left'
        self.loc_lr = 'lower right'   # location of legend
        self.loc_ll = 'lower left'    # location of legend
        self.loc_cr = 'center right'  # location of legend
        self.loc_cl = 'center left'   # location of legend
        
        # module class name
        mod_class = module.__class__.__module__
        self.mod_class_name = mod_class.split('.')[0]
        
        # continuing with SSC plots
        if self.mod_class_name == 'modules':
            
            # full PySAM module
            self.mod = module.Plant
            
            # define an Output object to extract information from SSC
            Outputs = self.mod.PySAM_Outputs if module.run_loop else self.mod.Outputs
            
            # saving full time logs
            self.t_full = np.asarray(Outputs.time_hr)*u.hr
            self.full_slice = slice(0, len(self.t_full), 1)
            
            # setting operating modes, kept it way at the bottom because it's ugly
            self.set_operating_modes_list()
            
            # extracting outputs
            SolarOutputExtraction.set_ssc_outputs(self, Outputs)
        
        
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
        power_array_list = ['p_cycle', 'q_dot_rec_in', 'q_thermal', 'gen', 'q_dot_pc_su'] # list of array strings
        power_label_list = ['P_cycle (Electric)',
                            'Q_dot Incident (Thermal)',
                            'Q_dot to Salt (Thermal)',
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
        ax2.set_ylim(-0.05*self.e_tes_design.m, 0.7*self.e_tes_design.m)

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


class SolarDispatchPlots(DispatchPlots):
    """
    The Plots class is a part of the PostProcessing family of classes. It can 
    output results from Pyomo Dispatch models. 

    Note that the DispatchPlots class must be initialized before using. 
    """

    def __init__(self, module, fsl='x-small', loc='best', legend_offset=False,
                 lp=16, lps=12, fs=12, lw=2, x_shrink=0.85, x_legend=12):
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
        SolarPlots.__init__( self, module, fsl, loc, legend_offset, \
                 lp, lps, fs, lw, x_shrink, x_legend)
 

        # continuing with Pyomo Dispatch plots
        if self.mod_class_name == 'dispatch':
            
            # saving module to self
            self.dm = module
            
            # extract outputs from Dispatch model
            SolarOutputExtraction.set_pyomo_outputs(self)
            
            # slice of arrays
            self.full_slice = slice(0, len(self.t_full), 1)


    def plot_pyomo_energy(self, ax=None, title_label=None, plot_all_time=True, \
                          start_hr=0, end_hr=48, hide_x=False,  x_legend=1.2, \
                          y_legend_L=1.0, y_legend_R=1.0):
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

        #========================#
        #--- Creating Figure  ---#
        #========================#

        # if no axis object specified, create a figure and axis from it
        if ax is None:
            fig = plt.figure(figsize=[10, 5])
            ax = fig.gca()   # this is the power plot

        # plot energy arrays
        energy_array_list = ['s_array', 'ucsu_array', 'ursu_array'] # list of array strings
        energy_label_list = ['TES Reserve Quantity',
                             'Cycle Startup Energy Inventory',
                             'Receiver Startup Energy Inventory'] # list of labels for each array string to extract from Outputs
        energy_wts = np.linspace(4, 2, 3).tolist()
        energy_ylabel = 'Energy \n(MWh)'

        ax  = self.plot_SSC_generic(ax, array_list=energy_array_list, \
                                    label_list=energy_label_list, \
                                    y_label=energy_ylabel, \
                                    lw_list=energy_wts, \
                                    title_label=title_label, \
                                    plot_all_time=plot_all_time, \
                                    start_hr=start_hr, end_hr=end_hr, hide_x=hide_x)
            
        #========================#
        #--- Extra arrays  ---#
        #========================#
        
        # retrieving arrays
        s_array = self.s_array.m    
        ucsu_array = self.ucsu_array.m
        ursu_array = self.ursu_array.m
        
        # vertical energy line at midpoint
        energy_vert = np.linspace( np.min([s_array, ucsu_array, ursu_array]), 
                                   np.max([s_array, ucsu_array, ursu_array])*1.1, 
                                   self.T )
            
        # Line marking the midpoint line
        ax.plot(self.time_midway, energy_vert, 'k--', linewidth=self.lw)


        #========================#
        #---- Setting Labels ----#
        #========================#
        
        # customizing legend(s)
        if self.legend_offset:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl, bbox_to_anchor=(x_legend, y_legend_L)) # plot legend for Energy arrays
        else:
            ax.legend( loc=self.loc_ul, fontsize=self.fsl) # plot legend for Energy arrays

        
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
                            'xrsu_array', 'wdot_s_array', 'wdot_p_array'] # list of array strings
        power_label_list = ['Cycle Out (E)', 'Cycle In (T)',
                            'Receiver Out (T)', 'Receiver Startup (T)',
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


    def plot_pyomo_solar_bin(self, ax=None, title_label=None, plot_all_time=True, \
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
        recbin_array_list = ['yr_array',   'yrhsp_array', 'yrsb_array',
                             'yrsd_array', 'yrsu_array',  'yrsup_array'] # list of array strings
        recbin_label_list = ['Is Receiver On?', 'Is Receiver HSU Pen?', 
                             'Is Receiver SB?', 'Is Receiver SD?', 
                             'Is Receiver SU?', 'Is Receiver CSU Pen?'] # list of labels for each array string to extract from Outputs
        recbin_wts = np.linspace(10, 1.5, 6).tolist()
        recbin_ylabel = 'Receiver \nBinary \nVariables'

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