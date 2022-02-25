#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 23:01:32 2022

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
from matplotlib.gridspec import GridSpec
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# =============================================================================
# json file names
# =============================================================================

jsons = ['model1_Hamilton_560_tariffx1',
         'model1_Hamilton_560_tariffx1_3',
         'model1_Hamilton_560_tariffx1_5',
         'model1_Hamilton_560_tariffx2'
         ]

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper")

# json name parameters
start_name     = 'failureModes'
PySAM_name     = 'PySAM' 
add_extra_Name = True
extra_name     = '2022_02' 
dispatch       = True 
sscH           = 24   
pyoH           = 48  
TES_min        = 0   
TES_max        = 7  
PC_min         = 450  
PC_max         = 1150  

# selecting coefficient array
coeffs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  

# key for coefficients
coeff_list = ['ec_p',     # proft term
              'ec_pcsu',  # PC cold SU term
              'ec_pchsu', # PC hot SU term 
              'ec_pcsd',  # PC SD term
              'ec_pcr',   # PC ramping term
              'ec_pcer',  # PC excess ramping term
              'ec_nsu',   # LFR cold SU termq
              'ec_nhsu',  # LFR hot SU term
              'ec_nsd',   # LFR SD term
              'ec_pcg',   # PC power generating term
              'ec_pcsb',  # PC SB term
              'ec_nt'     # LFR thermal power generating term
              ]

# coefficient strings
extr_str = ''
for n,coeff in enumerate(coeff_list):
    extr_str += "_{0}{1}".format( coeff.split('_')[1], str(coeffs[n]))

# generate name of file
def updated_json( json ):
    filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__{11}.nuctes'.format(
                    start_name,
                    PySAM_name,
                    json,
                    extra_name if add_extra_Name else '',
                    dispatch, sscH, pyoH,
                    TES_min, TES_max,
                    PC_min, PC_max, extr_str )
    return filename

# =============================================================================
# Loop through jsons
# =============================================================================

# ========== Figure ==========
fig = plt.figure(figsize=(14,10))
gs  = GridSpec(3,4, figure=fig)

plotrange = range(4)

for i in plotrange:
    ax1  = fig.add_subplot(gs[0:2,i])

    filename = updated_json( jsons[i] )
    
    # generate full path to file
    NTPath = os.path.join(output_dir, filename)
    
    # pickling
    with open(NTPath, 'rb') as f:
        Storage = pickle.load(f)
        
    # storage dictionary
    slicing = slice(5,15,1)
    
    tshours  = Storage['tshours'] 
    p_cycle  = Storage['p_cycle'][slicing] 
    ppa      = Storage['ppa'][:,slicing] 
    
    # =============================================================================
    # colormap
    # =============================================================================
    
    cmap=cm.seismic
    # ========== Arrays ==========
    array = copy.deepcopy( ppa )
    array[array==-1] = np.max(array)
    # array = np.array([array[n]/ppa.max(axis=1)[n] for n in range( len( ppa.max(axis=1) ) ) ])
    array /= array[0,-1]
    
    min_arr = 1 - array.min()
    max_arr = array.max() - 1
    max_diff = np.max([min_arr, max_arr])
    mean_label = "PPA Price Relative to Reference \n P=450MWe , TES=0"
    
    
    asp_df = 1.5
    im1 = ax1.imshow(array.T, origin='upper', cmap=cmap, aspect=asp_df, \
                     vmin = 1 - max_diff, vmax = 1 + max_diff )
    contours = plt.contour(array.T, levels=[0.95, 0.97, 0.99, 1.0], colors='black')
    ax1.clabel(contours, inline=1, fontsize=10)
    
    # ========== Text ==========
    xmin,xmax,ymax,ymin = im1.get_extent() #note the order
    x_series = np.linspace(xmin+0.5, xmax-0.5, len(tshours))
    y_series = np.linspace(ymin+0.5, ymax-0.5, len(p_cycle))
    
    # ========== labels ==========
    # setting axis labels
    ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
    
    # creating colorbar for the 2D heatmap with label
    # if i == plotrange[-1]:
    #     cb1 = fig.colorbar(im1, ax=ax1, pad=0.01)
    #     cb1.set_label(mean_label, labelpad= 8, fontweight = 'bold')
    
    # setting tick marks for x and y axes
    ax1.set_xticks(range(len(tshours)))
    ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )
    
    if i == 0:
        ax1.set_yticks(range(len(p_cycle)))
        ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle] )
        ax1.set_ylabel('Power Cycle Output\n(MWe)', fontweight='bold')
    else:
        ax1.set_yticklabels( '')


    # =============================================================================
    # 
    # =============================================================================
    # create object for NuclearTES for some setup steps
    nuctes = NuclearTES.NuclearTES( json_name= jsons[i], is_dispatch=False )
    default_tariff = nuctes.df_array
    
    
    # =============================================================================
    # Default Tariff Rate Plot
    # =============================================================================
    days_per_week = 7
    hrs_per_day   = 24
    weeks_till_summer = 26


    # indexing the winter default tariff schedule
    winter_end_slice = int(days_per_week * hrs_per_day)
    winter_slice = slice(0, winter_end_slice, 1)

    # indexing the summer default tariff schedule
    summer_start_slice = int( days_per_week * hrs_per_day * weeks_till_summer)
    summer_end_slice   = int( summer_start_slice + days_per_week * hrs_per_day)
    summer_slice = slice(summer_start_slice, summer_end_slice, 1)

    # time in units of days for plotting
    p_time = (np.arange(0,8760)*u.hr ).to('d').m

    # x labels
    xlabel_loc = np.arange(0, winter_end_slice, 24)
    xlabels = ["M", "Tu", "W", "Th", "F", "Sa", "Su"]

    #====== Figure ======#
    ax  = fig.add_subplot(gs[2,i])

    ax.plot(p_time[winter_slice], default_tariff[winter_slice], linewidth= 3, label="Non-Summer")
    ax.plot(p_time[winter_slice], default_tariff[summer_slice], linewidth= 3, label="Summer")

    miny = np.min(default_tariff)
    maxy = np.max(default_tariff)

    ax.set_ylim([0.2, 3.5])
    ax.set_xticks(p_time[winter_slice][xlabel_loc])
    ax.set_xticklabels(xlabels)
    ax.grid(True)

    # ax.set_xlabel("Time (d)", fontweight='bold')
    if i == 0:
        ax.set_ylabel("SAM Generic Peak \nPrice Multiplier", fontweight='bold')
        ax.legend()