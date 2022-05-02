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
import matplotlib.patheffects as pe
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy import interpolate
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# =============================================================================
# json file names
# =============================================================================

jsons = ['model1_Hamilton_560_tariffx1',
         # 'model1_Hamilton_560_tariffx1_3',
         'model1_Hamilton_560_tariffx1_5',
         'model1_Hamilton_560_tariffx2'
         ]

titles = ['Tariff Peaks x1',
          # 'Tariff Peaks x1.3',
          'Tariff Peaks x1.5',
          'Tariff Peaks x2']

nPlots = len(titles)

cbar_label = "PPA Price Relative to Reference \n @ P=450MWe , TES=0"

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper")

# json name parameters
start_name     = 'paramSweep_varTurbineCost' # failureModes
PySAM_name     = 'PySAM' 
add_extra_Name = True
extra_name     = '2022_03'  #  2022_02
dispatch       = True 
sscH           = 24   
pyoH           = 48  
TES_min        = 0   
TES_max        = 8    # 7
PC_min         = 450  
PC_max         = 900  

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
# for n,coeff in enumerate(coeff_list):
#     extr_str += "_{0}{1}".format( coeff.split('_')[1], str(coeffs[n]))

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
gs  = GridSpec(3,nPlots, figure=fig)

plotrange = range(nPlots)

for i in plotrange:
    ax1  = fig.add_subplot(gs[0:2,i])

    filename = updated_json( jsons[i] )
    
    # generate full path to file
    NTPath = os.path.join(output_dir, filename)
    
    # pickling
    with open(NTPath, 'rb') as f:
        Storage = pickle.load(f)
        
    # storage dictionary
    # slicing = slice(5,15,1)
    slicing = slice(0,15,1)
    
    tshours  = Storage['tshours'] 
    p_cycle  = Storage['p_cycle'][slicing] 
    ppa      = Storage['ppa'][:,slicing] 
    
    # =============================================================================
    # sanity checks
    # =============================================================================
    # no TES efficiency penalty
    print("************* {0} ******************** \n No TES Penalty".format(jsons[i]))
    opt_TES, opt_Pref = np.where( ppa == np.min(ppa) )
    print("Optimum PPA = {0} \n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa) / ppa[0,-1])*100, 
                                                                         p_cycle[opt_Pref][0], 
                                                                         tshours[opt_TES][0] )  )
    
    # some TES efficiency penalty
    pen = 1.01
    ppa_1 = copy.deepcopy(ppa)
    ppa_1[1:,:] *= pen
    print("************* {0} ******************** \n TES Penalty = {1}".format( jsons[i], (pen-1)*100 ) )
    opt_TES_p1, opt_Pref_p1 = np.where( ppa_1 == np.min(ppa_1) )
    print("Optimum PPA = {0} \n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa_1) / ppa[0,-1])*100,  
                                                                         p_cycle[opt_Pref_p1][0], 
                                                                         tshours[opt_TES_p1][0] )  )
    
    # some TES efficiency penalty
    pen = 1.025
    ppa_1 = copy.deepcopy(ppa)
    ppa_1[1:,:] *= pen
    print("************* {0} ******************** \n TES Penalty = {1}".format( jsons[i], (pen-1)*100 ) )
    opt_TES_p1, opt_Pref_p1 = np.where( ppa_1 == np.min(ppa_1) )
    print("Optimum PPA = {0} \n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa_1) / ppa[0,-1])*100,  
                                                                         p_cycle[opt_Pref_p1][0], 
                                                                         tshours[opt_TES_p1][0] )  )
    
    # some TES efficiency penalty
    pen = 1.05
    ppa_1 = copy.deepcopy(ppa)
    ppa_1[1:,:] *= pen
    print("************* {0} ******************** \n TES Penalty = {1}".format( jsons[i], (pen-1)*100 ) )
    opt_TES_p1, opt_Pref_p1 = np.where( ppa_1 == np.min(ppa_1) )
    print("Optimum PPA = {0} \n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa_1) / ppa[0,-1])*100,  
                                                                         p_cycle[opt_Pref_p1][0], 
                                                                         tshours[opt_TES_p1][0] )  )

    
    # =============================================================================
    # colormap
    # =============================================================================
    cmap=cm.seismic
    threshold = 1.109 # max threshold to plot
    threshold_contour = 1.1
    
    # ========== Arrays ==========
    
    # using PPA price as metric
    array = copy.deepcopy( ppa ) 
    # any values equal to -1 set to max (bad values, easier to spot)
    array[array==-1] = np.max(array)
    # normalizing by the reference @ 450 MWe and tshours = 0
    array /= array[0,-1]
    # anything above max threshold set to same value
    array[array>threshold] = threshold
    
    # flip array, since p_ref values are technically backwards
    array = np.fliplr(array)
    
    # getting upper and lower bounds for colorbar
    min_arr = 1 - array.min()
    max_arr = array.max() - 1
    max_diff = np.max([min_arr, max_arr])
    
    # plotting heatmap
    asp_df = 1.5
    im1 = ax1.imshow(array.T, origin='lower', #extent=[tshours[0], tshours[-1], p_cycle[-1], p_cycle[0] ],
                    cmap=cmap, aspect=asp_df, 
                     vmin = 1 - max_diff, vmax = 1 + max_diff )
    ax1.set_ylim([-0.5, 9.5])
    
    # defining contours
    finterp = interpolate.interp2d(tshours, np.flipud(p_cycle), array.T, kind='cubic')
    tt = np.linspace(tshours[0], tshours[-1], 1000)
    pp = np.linspace(p_cycle[-1], p_cycle[0], 1000)
    ff = finterp(tt,pp)
    
    # if i>0:
    # contours = plt.contour(array.T, levels=[0.95, 0.97, 0.99, 1.0], colors='black')
    # contours = plt.contour(array.T, levels=[0.94, 0.96, 0.98, 1.0], colors='black')
    # contours = plt.contour(array.T, levels=[0.88, 0.94, 1.0], colors='black')
    c_levels = [0.91, 0.96, 0.98, 1.0, threshold_contour] if i > 0 else [threshold_contour]
    c_colors = 'k'
    contours = ax1.contour(ff, extent=[tshours[0], tshours[-1], 0, 9 ],
                                    levels=c_levels, 
                                    colors=c_colors)
    
    # labels for each contour level
    c_fmt = {}
    c_empty = {}
    for l,s in zip(contours.levels, c_levels):
        c_fmt[l] = ' {:.2f} '.format(s)
        c_empty[l] = '      '
        
    # Add empty labels, just want to get locations for text
    ax1.clabel(contours, fmt=c_empty, inline=True, fontsize=10)
    
    # sometimes, number of levels doesnt correspond to actual number of contours drawn
    if len(contours.labelTexts) != len(contours.levels):
        ind = 0 - len(contours.labelTexts)   # get actual number of contours
        c_levels_out = contours.levels[ind:] # actual contours are at the back of the list
        c_fmt_out    = { k:c_fmt[k] for k in c_levels_out } # new dict
    else:
        c_levels_out = contours.levels
        c_fmt_out    = c_fmt
    
    # draw white text with black border for rel ppa price
    for textLabels,l,s in zip(contours.labelTexts, c_levels_out, c_fmt_out):
        loc = textLabels.get_position()
        xloc = loc[0] - 0.65
        yloc = loc[1] - 0.1
        txt = ax1.text(xloc, yloc, c_fmt[l], size=13, fontweight='bold', color='w')
        txt.set_path_effects([pe.withStroke(linewidth=3, foreground='k')])
        plt.draw()
        
    # ========== Text ==========
    xmin,xmax,ymax,ymin = im1.get_extent() #note the order
    x_series = np.linspace(xmin+0.5, xmax-0.5, len(tshours))
    y_series = np.linspace(ymin+0.5, ymax-0.5, len(p_cycle))
    
    # ========== labels ==========
    # setting axis labels
    ax1.set_xlabel('tshours\n(hr)', fontsize=14, fontweight='bold')
    ax1.set_title(titles[i], fontweight='bold')

    # setting tick marks for x and y axes
    ax1.set_xticks(range(len(tshours)))
    ax1.set_xticklabels( ['{0}'.format(t) for t in tshours], fontsize=14 )
    ax1.tick_params(axis='x',labelsize=14)
    ax1.tick_params(axis='y',labelsize=14)
    
    if i == 0:
        ax1.set_yticks(range(len(p_cycle)))
        ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in np.flipud(p_cycle)] )
        ax1.set_ylabel('Power Cycle Output\n(MWe)', fontsize=14, fontweight='bold')
    else:
        ax1.set_yticks(range(len(p_cycle)))
        ax1.set_yticklabels('')


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

    ax.plot(p_time[winter_slice], default_tariff[winter_slice], linewidth= 3, label="Winter")
    ax.plot(p_time[winter_slice], default_tariff[summer_slice], linewidth= 3, label="Summer")

    miny = np.min(default_tariff)
    maxy = np.max(default_tariff)

    ax.set_ylim([0.2, 3.5])
    ax.set_xticks(p_time[winter_slice][xlabel_loc])
    ax.set_xticklabels(xlabels, fontsize=14)
    ax.tick_params(axis='x',labelsize=14)
    ax.tick_params(axis='y',labelsize=14)

    ax.grid(True)

    # ax.set_xlabel("Time (d)", fontweight='bold')
    if i == 0:
        ax.set_ylabel("SAM Generic Peak \nPrice Multiplier", fontsize=14, fontweight='bold')
        ax.legend()
    else:
        ax.yaxis.set_ticklabels([])
    
plt.tight_layout()

fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.875, 0.4, 0.025, 0.55])
fig.colorbar(im1, cax=cbar_ax)
cbar_ax.set_ylabel(cbar_label, labelpad=18, fontsize=14, fontweight='bold')
cbar_ax.tick_params(axis='y',labelsize=14)


fig_name = 'tariff_exaggeration_effect_on_ppa.pdf'
fig.savefig( os.path.join(output_dir, fig_name), dpi=300 )