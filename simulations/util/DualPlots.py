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
import pint
u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
rc('axes', linewidth=2)
rc('font', weight='bold', size=12)
from util.PostProcessing import OutputExtraction
from util.PostProcessing import Plots
from util.SolarPlots import SolarOutputExtraction
from util.SolarPlots import SolarPlots
from util.PostProcessing import DispatchPlots


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