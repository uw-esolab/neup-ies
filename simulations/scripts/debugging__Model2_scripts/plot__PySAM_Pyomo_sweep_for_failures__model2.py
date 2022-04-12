#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 16:55:39 2022

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

start_name     = 'failureModes' # testDefocus # failureModes # failureModesNoCSPSB
PySAM_name     = 'PySAM'  # PySAM  # ''
add_extra_Name = True
extra_name    = '2022_04'  # 2021_10  # 2021_11 # 2021_12
json_name     = 'model2_Hamilton_560_tariffx1'   # model2_CAISO_Hamilton  # model2_Hamilton_560_TwoPeaks_x1 # model2_Hamilton_560_tariffx1
dispatch      = True # True # False
sscH          = 24   # 12 # 24
pyoH          = 48   # 24 # 48
q_dot_CSP     = 500
TES_min       = 0    # 0  # 2
TES_max       = 12   # 14
PC_min        = 600  # 100 # 300 # 400 # 550 # 600 # 900
PC_max        = 850  # 500 # 850 # 1150

# generate name of file
filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__CSP_{11}.nuctes'.format(
                start_name,
                PySAM_name,
                json_name,
                extra_name if add_extra_Name else '',
                dispatch, sscH, pyoH,
                TES_min, TES_max,
                PC_min, PC_max, q_dot_CSP )


# generate full path to file
NTPath = os.path.join(output_dir, filename)

if os.path.exists( NTPath ):
    print("File found at {0}".format(NTPath) )
else:
    print("File not found! Possible matches: ")
    file_list       = os.listdir( output_dir )
    relevant_sname  = '{0}_{1}__'.format( start_name, PySAM_name )
    relevant_filenames = np.array([f for f in file_list if relevant_sname in f])
    
    splits = np.array([ s.split(  relevant_sname )[1] for s in relevant_filenames ])
    splits = np.array([ s.split('__')[0] for s in splits])
    
    print_names = relevant_filenames[ splits == json_name ].tolist()
    for name in print_names:
        print( name )
    
# =============================================================================
# extract data
# =============================================================================

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
allOff_log   = Storage['allOff_log']
exceptions   = Storage['exceptions'] 
op_modes_list     = Storage['op_modes_list'] 
pyomo_bad_log     = Storage['pyomo_bad_log']
pyomo_bad_idx_log = Storage['pyomo_bad_idx_log'] 
if 'revenue' in Storage.keys():
    revenue = Storage['revenue']

# setting figure title
full_title = "PySAM - {0} tariffs {1} Pyomo - {2:.0f}/{3:.0f}hr horizons ".format(
                    json, "with" if dispatch else "without", sscH, pyoH)


# =============================================================================
# Simulation Progression
# =============================================================================
cmap = cm.hot_r

# ========== Arrays ==========
# p_array = copy.deepcopy( iter_log[:,:-1] ) / 24.
# array = copy.deepcopy( iter_log[:,:-1] ) / 8760. * 100
p_array = copy.deepcopy( iter_log ) / 24.
array = copy.deepcopy( iter_log ) / 8760. * 100
mean_label = "% of simulation completed"

# ========== Figure ==========
fig = plt.figure(figsize=(8,8))
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
         if color == 'w':
             color = 'C3' if allOff_log[i,j] > 0 else color
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
# 
# =============================================================================
cmap = cm.Reds

# ========== Arrays ==========
p_array = copy.deepcopy( allOff_log ) 
array = copy.deepcopy( allOff_log ) 
mean_label = "All OFF"

# ========== Figure ==========
fig = plt.figure(figsize=(8,8))
ax1  = fig.add_subplot(111)
fig.suptitle(full_title, fontweight='bold')

asp_df = 0.7
im1 = ax1.imshow(array.T, origin='upper', cmap=cmap, aspect=asp_df,
                 vmax = array.max() ) #, vmin = -1*array.max() )

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

# =============================================================================
# Defocus Fraction
# =============================================================================
# cmap = cm.hot_r

# # ========== Arrays ==========
# p_array = copy.deepcopy( defocus[:,:,0] ) 
# array = copy.deepcopy( defocus[:,:,0] ) 
# mean_label = "defocus fraction"

# # ========== Figure ==========
# fig = plt.figure(figsize=(8,14))
# ax1  = fig.add_subplot(111)
# fig.suptitle(full_title, fontweight='bold')

# asp_df = 0.7
# im1 = ax1.imshow(array.T, origin='upper', cmap=cmap, aspect=asp_df)

# # ========== Text ==========
# xmin,xmax,ymax,ymin = im1.get_extent() #note the order
# x_series = np.linspace(xmin+0.5, xmax-0.5, len(tshours))
# y_series = np.linspace(ymin+0.5, ymax-0.5, len(p_cycle))

# x, y = np.meshgrid(x_series, y_series)
# for i, x_val in enumerate(x_series):
#     for j, y_val in enumerate(y_series):
#          condition = fail_log[i,j]
#          color = 'y' if condition == 1 else 'g' if condition == 2 else 'm' if condition == 4 else 'w' 
#          color = 'c' if pyomo_bad_log[i,j] else color
#          text = ax1.text(x_val, y_val, "{0:.2f}".format(p_array[i,j]), color=color, fontsize=10, va='center', ha='center')
#          text.set_path_effects([PathEffects.Stroke(linewidth=3, foreground='black'), PathEffects.Normal()])
#          del text

# # ========== labels ==========
# # setting axis labels
# ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
# ax1.set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')

# # creating colorbar for the 2D heatmap with label
# cb1 = fig.colorbar(im1, ax=ax1, pad=0.01)
# cb1.set_label(mean_label, labelpad= 8, fontweight = 'bold')

# # setting tick marks for x and y axes
# ax1.set_xticks(range(len(tshours)))
# ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )

# ax1.set_yticks(range(len(p_cycle)))
# ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle] )

# =============================================================================
# 
# =============================================================================
# if 'revenue' in Storage.keys():
#     cmap=cm.Greens
#     # ========== Arrays ==========
#     p_array = copy.deepcopy( revenue[:,2:] ) /1e9
#     array = copy.deepcopy( revenue[:,2:] ) /1e9
#     mean_label = "Simple Revenue Metric \n($B)"
    
#     # ========== Figure ==========
#     fig = plt.figure(figsize=(8,14))
#     ax1  = fig.add_subplot(111)
#     fig.suptitle("Simple Revenue Metric (in B$)", fontweight='bold')
    
#     asp_df = 1
#     im1 = ax1.imshow(array.T, origin='upper', cmap=cmap, aspect=asp_df)
    
#     # ========== Text ==========
#     xmin,xmax,ymax,ymin = im1.get_extent() #note the order
#     x_series = np.linspace(xmin+0.5, xmax-0.5, len(tshours))
#     y_series = np.linspace(ymin+0.5, ymax-0.5, len(p_cycle)-2)
    
#     # x, y = np.meshgrid(x_series, y_series)
#     # for i, x_val in enumerate(x_series):
#     #     for j, y_val in enumerate(y_series):
#     #          condition = fail_log[i,j]
#     #          color = 'w'
#     #          text = ax1.text(x_val, y_val, "{0:.2f}".format(p_array[i,j]), color=color, fontsize=10, va='center', ha='center')
#     #          text.set_path_effects([PathEffects.Stroke(linewidth=3, foreground='black'), PathEffects.Normal()])
#     #          del text
    
#     # ========== labels ==========
#     # setting axis labels
#     ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
#     ax1.set_ylabel('Power Cycle Output\n(MWe)', fontweight='bold')
    
#     # creating colorbar for the 2D heatmap with label
#     cb1 = fig.colorbar(im1, ax=ax1, pad=0.01)
#     cb1.set_label(mean_label, labelpad= 8, fontweight = 'bold')
    
#     # setting tick marks for x and y axes
#     ax1.set_xticks(range(len(tshours)))
#     ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )
    
#     ax1.set_yticks(range(len(p_cycle)-2))
#     ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle[2:]] )
