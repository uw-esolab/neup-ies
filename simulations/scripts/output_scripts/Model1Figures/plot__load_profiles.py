#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:06:55 2022

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
from util.NuclearTESLoadProfiles import NuclearTESLoadProfiles
import matplotlib.pyplot as plt
import os, pint, time, copy, pickle
import numpy as np
import pyomo.environ as pe
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
#  Simulation Parameters
# =============================================================================

# modifying inputs
json = "model1_Hamilton_560_tariffx1"   # model1_CAISO_Hamilton # model1_Hamilton_560_tariffx2 # model1_Hamilton_560_tariffx1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 700  # (MW)
tshours = 2    # (hr)

# =============================================================================
#  Extracting Information
# =============================================================================

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper" )
filename = 'loadProfile_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_{4}__PC_{5}.nuctes'.format(
                json, dispatch, sscH, pyoH, tshours, Pref )
NTPath = os.path.join(output_dir, filename)

# creating plotting util
Plots = NuclearTESLoadProfiles( NTPath )
Storage = Plots.Storage


# =============================================================================
# create figure
# =============================================================================
gen_data_array   = Plots.gen_dict
tch_data_array   = Plots.e_ch_tes_dict

    
# create full figure

gen_fig = plt.figure(figsize=(18,7))
gen_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
             fontsize=16, fontweight='bold')

tch_fig = plt.figure(figsize=(18,7))
tch_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
             fontsize=16, fontweight='bold')

count = 0
for dayOfWeek in ["WeekDay", "WeekEnd"]:
    for season in ["Winter", "Summer"]:
    
        # defining subplot of full figure (2x2)
        subplot = "22" + str(count+1)
        
        # creating current axis
        gen_ax__ = gen_fig.add_subplot( eval(subplot) )
        tch_ax__ = tch_fig.add_subplot( eval(subplot) )


        # booleans for tariff overlay
        is_winter  = season is "Winter"
        is_weekday = dayOfWeek is "WeekDay"
        
        # create violin plot and tariffs on twined axis
        plotting_gen_array = gen_data_array[season][dayOfWeek].tolist()
        plotting_tch_array = tch_data_array[season][dayOfWeek].tolist()
        
        Plots.create_violin_plot(    gen_ax__, plotting_gen_array, v_color='C9' )
        Plots.create_violin_plot(    tch_ax__, plotting_tch_array, v_color='C4' )
        
        if 'tariffx' in json:
            Plots.create_tariff_overlay( gen_ax__, is_winter, is_weekday, s_color='C1')
            Plots.create_tariff_overlay( tch_ax__, is_winter, is_weekday, s_color='C1')
        
        # some basic plotting labels
        # common_plot_things( gen_ax__ )
        # common_plot_things( tch_ax__, is_gen=False )
        
        if count == 0 or count == 1:
            gen_ax__.set_title(season, fontsize=16, fontweight='bold')
            tch_ax__.set_title(season, fontsize=16, fontweight='bold')
        if count == 0 or count == 2:
            gen_ax__.set_ylabel(dayOfWeek + '\n$\\regular_{Power\ Generated \ (MWe)}}$', fontsize=16, fontweight='bold')
            tch_ax__.set_ylabel(dayOfWeek + '\n$\\regular_{Power\ Generated \ (MWe)}}$', fontsize=16, fontweight='bold')
    
        count += 1
        
gen_fig.tight_layout()
tch_fig.tight_layout()


# else:
#     gen_fig = plt.figure(figsize=(18,5))
#     gen_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
#                  fontsize=16, fontweight='bold')
    
#     tch_fig = plt.figure(figsize=(18,5))
#     tch_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
#                  fontsize=16, fontweight='bold')
    
#     gen_ax__ = gen_fig.add_subplot( 111 )
#     tch_ax__ = tch_fig.add_subplot( 111 )
    
#     create_violin_plot(    gen_ax__, gen_data.tolist(), v_color='C9' )
#     # create_tariff_overlay( gen_ax__, is_winter, is_weekday, s_color='C1')
    
#     create_violin_plot(    tch_ax__, etes_data.tolist(), v_color='C4' )
#     # create_tariff_overlay( tch_ax__, is_winter, is_weekday, s_color='C1')
    