#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:55:09 2022

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import modules.DualPlantTES as DualPlantTES
from util.FileMethods import FileMethods
from util.NuclearTESLoadProfiles import NuclearTESLoadProfiles
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os, pint, time, copy, pickle
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
#  Simulation Parameters
# =============================================================================

case = 9 #1,3,5,7,9

# modifying inputs
json = "model2_CAISO_Hamilton_mod"   # model1_CAISO_Hamilton # model1_Hamilton_560_tariffx2 # model1_Hamilton_560_tariffx1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)

    
if case==0:
    Pref    = 1100 # (MW)
    tshours = 1    # (hr)
elif case==1:
    Pref = 1300
    tshours = 5
elif case==2:
    Pref = 650
    tshours = 1
elif case==3:
    Pref = 750
    tshours = 6
elif case==4:
    Pref = 300
    tshours = 4
elif case==5:
    Pref = 375
    tshours = 6
elif case==6:
    Pref = 225
    tshours=4
elif case==7:
    Pref = 225
    tshours=10
elif case==8:
    Pref = 225
    tshours=4
elif case==9:
    Pref = 175
    tshours=12
# =============================================================================
#  Extracting Information
# =============================================================================

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "" )

# generating filename 
filename = 'loadProfile_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_{4}__PC_{5}.nuctes'.format(
                json, dispatch, sscH, pyoH, tshours, Pref )

# full filepath
NTPath = os.path.join(output_dir, filename)

# creating plotting util
Plots = NuclearTESLoadProfiles( NTPath )

# =============================================================================
#  Creating Figure
# =============================================================================

gen_data_array   = Plots.gen_dict
tch_data_array   = Plots.e_ch_tes_dict

title = "CAISO Market"
    
# create full figure
full_fig = plt.figure(figsize=(18,10))
title = "{0} \n P={1} MWe, TES={2} hrs".format(title,Pref,tshours)
full_fig.suptitle(title, fontsize=18, fontweight='bold')

gsQ  = gridspec.GridSpec(3, 1, height_ratios=[1,2,2], figure=full_fig)

gsQ0 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gsQ[0])
gsQ1 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gsQ[1])
gsQ2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gsQ[2])

weekday = full_fig.add_subplot(gsQ[1], frameon=False)
weekday.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

weekend = full_fig.add_subplot(gsQ[2], frameon=False)
weekend.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

gen_color = 'mediumseagreen'
tch_color = 'blueviolet'
s_color   = 'darkgreen'

for j, season in enumerate(["Winter", "Summer"]):
    
    is_winter  = season is "Winter"
    
    price_ax__ = full_fig.add_subplot( gsQ0[j], facecolor='white')
    
    Plots.create_tariff_overlay( price_ax__, is_winter, s_color=s_color)
    price_ax__.set_title(season, fontsize=18, fontweight='bold')
    
    for i, dayOfWeek in enumerate(["WeekDay", "WeekEnd"]):
    
        
        plot_index = (i,j)
        # creating current axis
        if i == 0:
            gen_ax__ = full_fig.add_subplot( gsQ1[0,j], facecolor='papayawhip')
            tch_ax__ = full_fig.add_subplot( gsQ1[1,j], facecolor='papayawhip')
        else:
            gen_ax__ = full_fig.add_subplot( gsQ2[0,j], facecolor='white')
            tch_ax__ = full_fig.add_subplot( gsQ2[1,j], facecolor='white')

        # booleans for tariff overlay
        
        is_weekday = dayOfWeek is "WeekDay"
        
        # create violin plot and tariffs on twined axis
        plotting_gen_array = gen_data_array[season][dayOfWeek].tolist()
        plotting_tch_array = tch_data_array[season][dayOfWeek].tolist()
        


        hide_x = False if i == 1 else True
        Plots.create_violin_plot(  gen_ax__, plotting_gen_array, gen_data_array,v_color=gen_color, hide_x=True )
        Plots.create_violin_plot(  tch_ax__, plotting_tch_array, tch_data_array, v_color=tch_color, hide_x=hide_x )

        gen_ax__.xaxis.set_ticklabels([])
        
        if i == 0:
            weekday.set_ylabel("Weekdays", labelpad=80, fontsize=18, fontweight='bold')
        else:
            weekend.set_ylabel("Weekends", labelpad=80, fontsize=18, fontweight='bold')
        
        if j == 0:
            gen_ax__.set_ylabel("Power \nGenerated \n(MWe)", labelpad=16, fontsize=12, fontweight='bold')
            tch_ax__.set_ylabel("Tank \nCharge \n(MWh)", fontsize=12, fontweight='bold')

        
full_fig.tight_layout()


