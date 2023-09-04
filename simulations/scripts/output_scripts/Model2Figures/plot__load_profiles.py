#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:06:55 2022

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

case =11

# modifying inputs
if case%2 == 0:
    json = "model2_Hamilton_560_tariffx1_mod"
else:
    raise ValueError("Use other plotting script")

dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)


    
if case==0:
    Pref    = 1160 # (MW)
    tshours = 0.22    # (hr)
elif case==1:
    Pref = 1940
    tshours = 4.87
elif case==2:
    Pref = 710
    tshours = 0.37
elif case==3:
    Pref = 1040
    tshours = 5.62
elif case==4:
    Pref = 378
    tshours = 0.69
elif case==5:
    Pref = 363.7
    tshours = 8.61
elif case==6:
    Pref = 307
    tshours=0.85
elif case==7:
    Pref = 229.5
    tshours=11.32
elif case==8:
    Pref = 269.5
    tshours=0.96
elif case==9:
    Pref = 154.7
    tshours=14.86
elif case==10:
    Pref=260
    tshours=1
elif case==11:
    Pref=140
    tshours=16

    
# =============================================================================
#  Extracting Information
# =============================================================================

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "" )

# generating filename 
filename = 'loadProfile_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_{4}__PC_{5}_{6}.nuctes'.format(
                json, dispatch, sscH, pyoH, tshours, Pref,case )

# full filepath
NTPath = os.path.join(output_dir, filename)

# creating plotting util
Plots = NuclearTESLoadProfiles( NTPath )

# =============================================================================
#  Creating Figure
# =============================================================================

gen_data_array   = Plots.gen_dict
tch_data_array   = Plots.e_ch_tes_dict

title = "SAM Generic Peak Market {0}".format("x1" if "x1" in json else "x2")
    
# create full figure
full_fig = plt.figure(figsize=(18,10))
title = "{0} \n P={1} MWe, TES={2} hrs".format(title,Pref,tshours)
full_fig.suptitle(title, fontsize=18, fontweight='bold')

gsQ  = gridspec.GridSpec(2,1, figure=full_fig)

gsQ1 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gsQ[0])
gsQ2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gsQ[1])

weekday = full_fig.add_subplot(gsQ[0], frameon=False)
weekday.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

weekend = full_fig.add_subplot(gsQ[1], frameon=False)
weekend.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

gen_color = 'mediumseagreen'
tch_color = 'blueviolet'
s_color   = 'C1' if 'tariffx' in json else 'C0'

count = 0
for i, dayOfWeek in enumerate(["WeekDay", "WeekEnd"]):
    for j, season in enumerate(["Winter", "Summer"]):
        
        plot_index = (i,j)
        # creating current axis
        if i == 0:
            gen_ax__ = full_fig.add_subplot( gsQ1[0,j], facecolor='white')
            tch_ax__ = full_fig.add_subplot( gsQ1[1,j], facecolor='white')
        else:
            gen_ax__ = full_fig.add_subplot( gsQ2[0,j], facecolor='papayawhip')
            tch_ax__ = full_fig.add_subplot( gsQ2[1,j], facecolor='papayawhip')

        # booleans for tariff overlay
        is_winter  = season is "Winter"
        is_weekday = dayOfWeek is "WeekDay"
        
        # create violin plot and tariffs on twined axis
        plotting_gen_array = gen_data_array[season][dayOfWeek].tolist()
        plotting_tch_array = tch_data_array[season][dayOfWeek].tolist()
        
        Plots.create_tariff_overlay( gen_ax__, is_winter, is_weekday, s_color=s_color)
        Plots.create_tariff_overlay( tch_ax__, is_winter, is_weekday, s_color=s_color)

        hide_x = False if i == 1 else True
        Plots.create_violin_plot(  gen_ax__, plotting_gen_array, gen_data_array, v_color=gen_color, hide_x=True )
        Plots.create_violin_plot(  tch_ax__, plotting_tch_array, tch_data_array, v_color=tch_color, hide_x=hide_x )

        gen_ax__.xaxis.set_ticklabels([])
        
        if i == 0:
            gen_ax__.set_title(season, fontsize=18, fontweight='bold')
            weekday.set_ylabel("Weekdays", labelpad=80, fontsize=18, fontweight='bold')
        else:
            weekend.set_ylabel("Weekends", labelpad=80, fontsize=18, fontweight='bold')
        
        if j == 0:
            gen_ax__.set_ylabel("Power \nGenerated \n(MWe)", labelpad=16, fontsize=14, fontweight='bold')
            tch_ax__.set_ylabel("Tank \nCharge \n(MWh)", fontsize=14, fontweight='bold')

        count += 1
        
full_fig.tight_layout()


