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

case = 7

if case%2 == 0:
    jsons = ['model2_Hamilton_560_tariffx1_mod']
else:
    jsons = ['model2_CAISO_Hamilton_mod']
    
titles = ['Price Multipliers']



# locating output directory
output_dir = os.path.join( FileMethods.output_dir,"")

# json name parameters
start_name     = 'paramSweep_varTurbineCost' # failureModes
PySAM_name     = 'PySAM' 
add_extra_Name = True
extra_name     = '2022_05'  #  2022_05
dispatch       = True 
sscH           = 24   
pyoH           = 48  
solar_json = "__115"

if case <=1:
    extr_str = 'zero'
elif case <=3:
    extr_str = 'sweep'
elif case <=5:
    extr_str = 'small1'
elif case <=7:
    extr_str = 'small2'
elif case <=9:
    extr_str = 'large'
elif case <=11:
    extr_str = 'micro'


"""
elif case <=9:
    extr_str = 'verym'
elif case<=11:
    extr_str = 'nano'
else:
    extr_str = 'nearz'
"""

if extr_str == 'sweep':
    PC_min= 550
    PC_max=1000
    TES_min        = 0   
    TES_max        = 10   
    q_dot_nuclear_des= 950


elif extr_str == 'large':
    PC_min         = 1000  
    PC_max         = 1900  
    TES_min        = 0   
    TES_max        = 10   
    q_dot_nuclear_des = 1900

elif extr_str == 'small1':
    PC_min         = 150  
    PC_max         = 600  
    TES_min        = 0   
    TES_max        = 16   
    q_dot_nuclear_des=250

elif extr_str == 'small2':
    PC_min         = 125  
    PC_max         = 350 
    TES_min        = 0   
    TES_max        = 16  
    q_dot_nuclear_des=100

    
elif extr_str == 'zero':
    PC_min         = 115
    PC_max         = 220  
    TES_min        = 0   
    TES_max        = 16 

elif extr_str == 'micro':
    PC_min = 120
    PC_max = 300
    TES_min = 0
    TES_max = 16
    q_dot_nuclear_des=20
    
if solar_json=="__1":
    PC_min = int(q_dot_nuclear_des*450/950)
    if case % 2 ==0:
        PC_max=PC_min
        TES_max=0
    else:
        PC_max=2*PC_min
        TES_max=10
    

cbar_label = "PPA Price (c/kWe)"

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

# for n,coeff in enumerate(coeff_list):
#     extr_str += "_{0}{1}".format( coeff.split('_')[1], str(coeffs[n]))

# generate name of file
def updated_json( json ):
    filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__{11}{12}.nuctes'.format(
                    start_name,
                    PySAM_name,
                    json,
                    extra_name if add_extra_Name else '',
                    dispatch, sscH, pyoH,
                    TES_min, TES_max,
                    PC_min, PC_max, extr_str,solar_json )
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


if case==6:
    #issues with zero thermal storage
    ppa=np.delete(ppa,0,axis=0)
    ppa=np.delete(ppa,-1,axis=1)
    tshours=tshours[1:]
    p_cycle=p_cycle[:-1]

if case==10:
    for j in range(3):
        ppa=np.delete(ppa,-1,axis=1)
    p_cycle=p_cycle[:-3]



# ppa[1:,:] *= 1.05

# =============================================================================
# sanity checks
# =============================================================================
# no TES efficiency penalty
print("************* {0} ******************** \n No TES Penalty".format(jsons[0]))
opt_TES, opt_Pref = np.where( ppa == np.min(ppa) )


print("Optimum PPA = {3}, reduction = {0} %\n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa) / ppa[0,-1])*100, 
                                                                     p_cycle[opt_Pref][0], 
                                                                     tshours[opt_TES][0],
                                                                     np.min(ppa))  )
print('Cap fac = {0}'.format(Storage['cap_fac'][opt_TES[0]][opt_Pref[0]]))

# some TES efficiency penalty
pen = 1.01
ppa_1 = copy.deepcopy(ppa)
ppa_1[1:,:] *= pen
print("************* {0} ******************** \n TES Penalty = {1}".format( jsons[0], (pen-1)*100 ) )
opt_TES_p1, opt_Pref_p1 = np.where( ppa_1 == np.min(ppa_1) )
print("Optimum PPA = {3}, reduction = {0} %\n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa_1) / ppa[0,-1])*100,  
                                                                     p_cycle[opt_Pref_p1][0], 
                                                                     tshours[opt_TES_p1][0],
                                                                     np.min(ppa_1))  )

# some TES efficiency penalty
pen = 1.025
ppa_1 = copy.deepcopy(ppa)
ppa_1[1:,:] *= pen
print("************* {0} ******************** \n TES Penalty = {1}".format( jsons[0], (pen-1)*100 ) )
opt_TES_p1, opt_Pref_p1 = np.where( ppa_1 == np.min(ppa_1) )
print("Optimum PPA = {3}, reduction = {0} %\n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa_1) / ppa[0,-1])*100,  
                                                                     p_cycle[opt_Pref_p1][0], 
                                                                     tshours[opt_TES_p1][0],
                                                                     np.min(ppa_1))  )

# some TES efficiency penalty
pen = 1.11
ppa_1 = copy.deepcopy(ppa)
ppa_1[1:,:] *= pen
print("************* {0} ******************** \n TES Penalty = {1}".format( jsons[0], (pen-1)*100 ) )
opt_TES_p1, opt_Pref_p1 = np.where( ppa_1 == np.min(ppa_1) )
print("Optimum PPA = {3}, reduction = {0} %\n ****** at P = {1}MWe and TES = {2}hrs \n".format( (1 - np.min(ppa_1) / ppa[0,-1])*100,  
                                                                     p_cycle[opt_Pref_p1][0], 
                                                                     tshours[opt_TES_p1][0],
                                                                     np.min(ppa_1))  )

# =============================================================================
# colormap
# =============================================================================
cmap=cm.seismic
threshold = 1.12 # max threshold to plot
threshold_contour = np.round(1.1,2)

if case==7:
    ppa[3,3]=np.nan 

# ========== Arrays ==========

# using PPA price as metric
array = copy.deepcopy( ppa ) 
    
# any values equal to -1 set to max (bad values, easier to spot)
array[array==-1] = np.max(array)
# normalizing by the reference @ 450 MWe and tshours = 0
#array /= array[0,-1]
# anything above max threshold set to same value
#array[array>threshold] = threshold

# flip array, since p_ref values are technically backwards
array = np.fliplr(array)


# getting upper and lower bounds for colorbar
min_arr = 1 - array.min()
max_arr = array.max() - 1
max_diff = np.max([min_arr, max_arr])

# plotting heatmap
asp_df = 1
im1 = ax.imshow(array.T, origin='lower', #extent=[tshours[0], tshours[-1], p_cycle[-1], p_cycle[0] ],
                cmap=cmap, aspect=asp_df)#, 
                 #vmin = 1 - max_diff, vmax = 1 + max_diff )\

if extr_str == 'large' or extr_str == 'sweep':
    ax.set_ylim([-0.5, 9.5])
elif extr_str == 'small2':
    ax.set_ylim([-0.5,8.5])
    ax.set_xlim([-0.5,7.5])
else:
    ax.set_ylim([-0.5,6.5])
    ax.set_xlim([-0.5,8.5])
    
# defining contours
finterp = interpolate.interp2d(tshours, np.flipud(p_cycle), array.T, kind='cubic')
tt = np.linspace(tshours[0], tshours[-1], 1000)
pp = np.linspace(p_cycle[-1], p_cycle[0], 1000)
ff = finterp(tt,pp)

# contours = plt.contour(array.T, levels=[0.88, 0.94, 1.0], colors='black')
c_levels = [0.85, 0.9, 0.95, 1.0, threshold_contour]
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

fig.subplots_adjust(right=0.8 if extr_str in ['sweep','large'] else 0.9 if extr_str == 'small2' else 0.75)
cbar_ax = fig.add_axes([0.8, 0.1, 0.035, 0.875])
fig.colorbar(im1, cax=cbar_ax)
cbar_ax.set_ylabel(cbar_label, labelpad=20, fontweight='bold', fontsize=16)
