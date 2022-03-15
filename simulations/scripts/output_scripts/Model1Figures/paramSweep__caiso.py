#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:30:12 2022

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
from scipy import interpolate
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=16)
u = pint.UnitRegistry()

# =============================================================================
# json file names
# =============================================================================

jsons = ['model1_CAISO_Hamilton']

titles = ['CAISO Price Multipliers']

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
PC_max         = 900  # 1150

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
fig = plt.figure(figsize=(11,8))
ax  = fig.add_subplot(111)

# get CAISO 
filename = updated_json( jsons[0] )
    
# generate full path to file
NTPath = os.path.join(output_dir, filename)
    
# pickling
with open(NTPath, 'rb') as f:
    Storage = pickle.load(f)
        
# storage dictionary
slicing = slice(0,15,1) # full slice 
# slicing = slice(5,15,1) # partial slice over Pref

# extract data
tshours  = Storage['tshours'] 
p_cycle  = Storage['p_cycle'][slicing] 
ppa      = Storage['ppa'][:,slicing] 

# =============================================================================
# sanity checks
# =============================================================================

opt_TES, opt_Pref = np.where( ppa == np.min(ppa) )
print("Optimum PPA = {0} \n ****** at P = {1}MWe and TES = {2}hrs".format( np.min(ppa), 
                                                                     p_cycle[opt_Pref][0], 
                                                                     tshours[opt_TES][0] )  )


# =============================================================================
# colormap
# =============================================================================
cmap=cm.seismic
threshold = 1.12 # max threshold to plot
threshold_contour = np.round(1.1,2)

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
asp_df = 1
im1 = ax.imshow(array.T, origin='lower', #extent=[tshours[0], tshours[-1], p_cycle[-1], p_cycle[0] ],
                cmap=cmap, aspect=asp_df, 
                 vmin = 1 - max_diff, vmax = 1 + max_diff )
ax.set_ylim([-0.5, 9.5])
    
# defining contours
finterp = interpolate.interp2d(tshours, np.flipud(p_cycle), array.T, kind='cubic')
tt = np.linspace(tshours[0], tshours[-1], 1000)
pp = np.linspace(p_cycle[-1], p_cycle[0], 1000)
ff = finterp(tt,pp)

# contours = plt.contour(array.T, levels=[0.88, 0.94, 1.0], colors='black')
c_levels = [0.91, 0.94, 0.97, 1.0, threshold_contour]
c_colors = 'k'
contours = plt.contour(ff, extent=[tshours[0], tshours[-1], 0, 9 ],
                                levels=c_levels, 
                                colors=c_colors)

# labels for each contour level
c_fmt = {}
for l,s in zip(contours.levels, c_levels):
    c_fmt[l] = ' {:.2f} '.format(s)
    
# Add contour labels
ax.clabel(contours, fmt=c_fmt, inline=1, fontsize=16)
    

# ========== labels ==========
# setting axis labels
ax.set_xlabel('tshours\n(hr)', fontweight='bold', fontsize=16)

# setting tick marks for x and y axes
ax.set_xticks(range(len(tshours)))
ax.set_xticklabels( ['{0}'.format(t) for t in tshours] )
    

ax.set_yticks(range(len(p_cycle)))
ax.set_yticklabels( ['{0:.0f}'.format(P) for P in np.flipud( p_cycle ) ] )
ax.set_ylabel('Power Cycle Output\n(MWe)', fontweight='bold', fontsize=16)


plt.tight_layout()

fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.8, 0.1, 0.035, 0.875])
fig.colorbar(im1, cax=cbar_ax)
cbar_ax.set_ylabel(cbar_label, labelpad=20, fontweight='bold', fontsize=16)
