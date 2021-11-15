#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 22:42:41 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
from util.FileMethods import FileMethods
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)


# =============================================================================
# Retrieving Data
# =============================================================================

# locating output directory
output_dir = FileMethods.output_dir
# filename   = 'testDefocus__model1_noMin__pyomo_0__horizon_24_48__TES_[0,14]__PC_[400,850].nuctes' 
# filename   = 'testDefocus__model1__pyomo_0__horizon_24_48__TES_[0,14]__PC_[400,850].nuctes' 
# filename   = 'testDefocus__model1__pyomo_1__horizon_24_48__TES_[0,14]__PC_[400,850].nuctes'
# filename   = 'testDefocus__model1_CAISO__pyomo_0__horizon_24_48__TES_[0,14]__PC_[400,850].nuctes' 
# filename   = 'testDefocus__model1_CAISO__pyomo_1__horizon_24_48__TES_[0,14]__PC_[400,850].nuctes' 
# filename   = 'testDefocus__model1_CAISO__pyomo_1__horizon_12_24__TES_[0,14]__PC_[400,850].nuctes' 

# filename   = 'failureModes__model1_2021_10__pyomo_0__horizon_24_48__TES_[0,14]__PC_[100,850].nuctes' 
# filename   = 'failureModes__model1_2021_11__pyomo_0__horizon_24_48__TES_[0,14]__PC_[100,850].nuctes' 

# filename   = 'failureModes_PySAM__model1_2021_11__pyomo_0__horizon_24_48__TES_[0,14]__PC_[100,850].nuctes' 
# filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_0__horizon_24_48__TES_[0,14]__PC_[100,850].nuctes'
# filename   = 'failureModes_PySAM__model1_2021_11__pyomo_1__horizon_24_48__TES_[0,14]__PC_[100,850].nuctes'
# filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_1__horizon_24_48__TES_[0,14]__PC_[100,850].nuctes'
# filename   = 'failureModes_PySAM__model1_2021_11__pyomo_1__horizon_12_24__TES_[0,14]__PC_[100,850].nuctes'
# filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_1__horizon_12_24__TES_[0,14]__PC_[100,850].nuctes'

# filename   = 'failureModes_PySAM__model1_2021_11__pyomo_1__horizon_24_48__TES_[0,14]__PC_[100,500].nuctes'
# filename   = 'failureModes_PySAM__model1_2021_11__pyomo_1__horizon_12_24__TES_[0,14]__PC_[100,500].nuctes'
# filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_1__horizon_24_48__TES_[0,14]__PC_[100,500].nuctes'
# filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_1__horizon_12_24__TES_[0,14]__PC_[100,500].nuctes'
# filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_1__horizon_24_48__TES_[0,14]__PC_[550,850].nuctes'
filename   = 'failureModes_PySAM__model1_CAISO_2021_11__pyomo_1__horizon_12_24__TES_[0,14]__PC_[550,850].nuctes'

NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'rb') as f:
    Storage = pickle.load(f)
    
# sotrage dictionary
time_elapsed = Storage['time_elapsed']
dispatch     = Storage['dispatch']  
tshours      = Storage['tshours'] 
p_cycle      = Storage['p_cycle'] 
iterator1    = Storage['iterator1']  
iterator2    = Storage['iterator2']  
sscH         = Storage['sscH']   
pyoH         = Storage['pyoH']  
json         = Storage['json'] 
run_loop     = Storage['run_loop']   
defocus      = Storage['defocus']  
gen          = Storage['gen'] 
qdot_nuc     = Storage['qdot_nuc']   
TES_CH       = Storage['TES_CH']  
TES_DC       = Storage['TES_DC'] 
iter_log     = Storage['iter_log']   
fail_log     = Storage['fail_log']  
exceptions   = Storage['exceptions'] 
op_modes_list     = Storage['op_modes_list'] 
pyomo_bad_log     = Storage['pyomo_bad_log']
pyomo_bad_idx_log = Storage['pyomo_bad_idx_log'] 

# setting figure title
full_title = "PySSC - {0} tariffs {1} Pyomo - {2:.0f}/{3:.0f}hr horizons ".format(
                    json, "with" if dispatch else "without", sscH, pyoH)


# =============================================================================
# Simulation Progression
# =============================================================================
cmap = cm.hot_r

# ========== Arrays ==========
p_array = copy.deepcopy( iter_log ) / 24.
array = copy.deepcopy( iter_log ) / 8760. * 100
mean_label = "% of simulation completed"

# ========== Figure ==========
fig = plt.figure(figsize=(8,14))
ax1  = fig.add_subplot(111)
fig.suptitle(full_title, fontweight='bold')

asp_df = 0.7
im1 = ax1.imshow(array.T, origin='upper', cmap=cmap, aspect=asp_df, vmin=0, vmax=100)

# ========== Text ==========
xmin,xmax,ymax,ymin = im1.get_extent() #note the order
x_series = np.linspace(xmin+0.5, xmax-0.5, len(tshours))
y_series = np.linspace(ymin+0.5, ymax-0.5, len(p_cycle))

x, y = np.meshgrid(x_series, y_series)
for i, x_val in enumerate(x_series):
    for j, y_val in enumerate(y_series):
         condition = fail_log[i,j]
         color = 'y' if condition == 1 else 'g' if condition == 2 else 'm' if condition == 4 else 'w' 
         color = 'c' if pyomo_bad_log[i,j] else color
         text = ax1.text(x_val, y_val, '{0:.1f}'.format(p_array[i,j]), color=color, fontsize=10, va='center', ha='center')
         text.set_path_effects([PathEffects.Stroke(linewidth=3, foreground='black'), PathEffects.Normal()])
         del text

# ========== labels ==========
# setting axis labels
ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
ax1.set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')

# creating colorbar for the 2D heatmap with label
cb1 = fig.colorbar(im1, ax=ax1, pad=0.01)
cb1.set_label(mean_label, labelpad= 8, fontweight = 'bold')

# setting tick marks for x and y axes
ax1.set_xticks(range(len(tshours)))
ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )

ax1.set_yticks(range(len(p_cycle)))
ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle] )


# =============================================================================
# Defocus Fraction
# =============================================================================
cmap = cm.hot_r

# ========== Arrays ==========
p_array = copy.deepcopy( defocus[:,:,0] ) 
array = copy.deepcopy( defocus[:,:,0] ) 
mean_label = "defocus fraction"

# ========== Figure ==========
fig = plt.figure(figsize=(8,14))
ax1  = fig.add_subplot(111)
fig.suptitle(full_title, fontweight='bold')

asp_df = 0.7
im1 = ax1.imshow(array.T, origin='upper', cmap=cmap, aspect=asp_df)

# ========== Text ==========
xmin,xmax,ymax,ymin = im1.get_extent() #note the order
x_series = np.linspace(xmin+0.5, xmax-0.5, len(tshours))
y_series = np.linspace(ymin+0.5, ymax-0.5, len(p_cycle))

x, y = np.meshgrid(x_series, y_series)
for i, x_val in enumerate(x_series):
    for j, y_val in enumerate(y_series):
         condition = fail_log[i,j]
         color = 'y' if condition == 1 else 'g' if condition == 2 else 'm' if condition == 4 else 'w' 
         color = 'c' if pyomo_bad_log[i,j] else color
         text = ax1.text(x_val, y_val, "{0:.2f}".format(p_array[i,j]), color=color, fontsize=10, va='center', ha='center')
         text.set_path_effects([PathEffects.Stroke(linewidth=3, foreground='black'), PathEffects.Normal()])
         del text

# ========== labels ==========
# setting axis labels
ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
ax1.set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')

# creating colorbar for the 2D heatmap with label
cb1 = fig.colorbar(im1, ax=ax1, pad=0.01)
cb1.set_label(mean_label, labelpad= 8, fontweight = 'bold')

# setting tick marks for x and y axes
ax1.set_xticks(range(len(tshours)))
ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )

ax1.set_yticks(range(len(p_cycle)))
ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle] )
