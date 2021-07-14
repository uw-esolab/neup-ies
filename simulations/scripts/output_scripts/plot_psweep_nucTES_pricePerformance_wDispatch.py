#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 09:58:26 2021

@author: gabrielsoto
"""

import os
import matplotlib.pyplot as plt
import pint
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()
import pickle
from util.FileMethods import FileMethods
import matplotlib.cm as cm

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# Retrieving Data
# =============================================================================

# locating output directory
output_dir = FileMethods.output_dir
filename   = 'pricePerfvsDispatch_sizingTESandCycle.nuctes' 
NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'rb') as f:
    Storage = pickle.load(f)

# sotrage dictionary
time_elapsed = Storage['time_elapsed']
dispatch     = Storage['dispatch']  
tshours      = Storage['tshours'] 
p_mult       = Storage['p_mult'] 
P_ref        = Storage['P_ref']  
iterator1    = Storage['iterator1']  
iterator2    = Storage['iterator2']  
sscH         = Storage['sscH']   
pyoH         = Storage['pyoH']   
annual_energy_array = Storage['annual_energy_array'] 
ppa_array           = Storage['ppa_array']
lcoe_nom_array      = Storage['lcoe_nom_array']

    
# =============================================================================
# Set Up Figures and Initialize Lists
# =============================================================================

# creating figure
fig = plt.figure(figsize=[18.5,6])
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# adjusting space between 
plt.subplots_adjust(wspace=0.4)

# names of dispatch scenarios
dispatch_scenarios  = ['Pyomo - 24hr Horizon', 'Pyomo - 48hr Horizon', 'No Pyomo - SSC Only ']  

# specific performance metrics to plot 
performance_metrics = [annual_energy_array, ppa_array, lcoe_nom_array]
N              = len(performance_metrics)
get_max        = [True, False, False] # whether each performance is optimized with max or min
axes           = [ax1,ax2,ax3] # list of axes
heat_cmap_list = ['Blues', 'Blues_r', 'Blues_r'] # colormaps used for each performance metric
metric_labels  = ['Annual Energy Production \n(TWh)',
                  'PPA Price \n(cents/kWh)',
                  'LCOE \n(cents/kWh)' ] # labels for each performance metric

# colormap for contours
cmap_contour = cm.RdBu

# some other arrays
P_array  = P_ref * p_mult #actual power cycle reference outputs
lp = 12 #labelpad


# ========================================   
# Looping!
# ========================================     

for n in range(N):
    
    # indeces for dispatch scenarios
    dispatch_series  = list( range(N) ) # declare this every loop
    
    #extracting specific performance metric
    performance_metric = performance_metrics[n]
    
    # get optimum value of performance metric 
    Opt = np.amax(performance_metric) if get_max[n] \
            else np.amin(performance_metric)

    # indeces for this maximum value - { D=Dispatch_Type, T=tshours, P=power_cycle }
    [Opt_D, Opt_T, Opt_P] = [np.where(performance_metric == Opt)[d][0] for d in range(N)]

    # plotting 2D heat map of the **Dispatch type** that leads to optimum value
    im = axes[n].imshow(performance_metric[Opt_D].T, origin='lower', 
                        cmap=heat_cmap_list[n])

    # setting axis labels
    axes[n].set_xlabel('tshours \n(hr)', fontweight='bold')
    if n == 0:
        axes[n].set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')

    # set title of plot
    title = axes[n].set_title('Optimum: ' + dispatch_scenarios[Opt_D], fontweight='bold')

    # creating colorbar for the 2D heatmap with label
    cb = fig.colorbar(im, ax=axes[n], fraction=0.06, pad=0.01)
    cb.set_label(metric_labels[n], labelpad= lp, fontweight = 'bold')

    # setting tick marks for x and y axes
    axes[n].set_xticks(range(len(tshours)))
    axes[n].set_xticklabels( ['{0}'.format(t) for t in tshours] )
    if n == 0:
        axes[n].set_yticks(range(len(P_array)))
        axes[n].set_yticklabels( ['{0:.2f}'.format(P) for P in P_array] )
    else:
        axes[n].axes.yaxis.set_visible(False)

    # getting extent of heatmap to help with contours
    xmin,xmax,ymin,ymax = im.get_extent()
    x_series = np.linspace(xmin, xmax, len(tshours))
    y_series = np.linspace(ymin, ymax, len(P_array))

    #=========================================================================
    ### First Contour: Difference between Optimum and Suboptimum #1
    rm = dispatch_series.pop(Opt_D) #remove optimal index
    
    # find difference between optimal Dispatch type and next available
    diff_1 = performance_metric[Opt_D] - performance_metric[dispatch_series[0]]
    diff_1 /= performance_metric[Opt_D]
    
    # find min and max values of difference
    diff1min = diff_1.min()
    diff1max = diff_1.max()
    
    # define contour levels
    if np.sign(diff1min) != np.sign(diff1max):
        # difference has only one sign, so optimal is strictly better
        levels_1 = [diff1min, 0.7*diff1min, 0.4*diff1min, 0, 0.4*diff1max, 0.9*diff1max, diff1max]
    else:
        # difference changes signs, so sometimes the Dispatch type behaves worse
        levels_1 = [diff1min, 0.75*diff1min, 0.3*diff1min, diff1max]

    # contour plot for Suboptimum #1
    CS1 = axes[n].contour(x_series, y_series, diff_1.T, levels_1, cmap=cmap_contour)

    # labels for each contour level
    fmt_1 = {}
    for l,s in zip(CS1.levels, levels_1):
        # plotting as percentages
        fmt_1[l] = ' {0:.1f}'.format(s*100) + ' % '
        
    # Add contour labels
    axes[n].clabel(CS1,  inline=True, fmt=fmt_1, fontsize=10)
    
    # creating label for entire contour plot
    CS1.collections[n].set_label('Pct wrt ' + dispatch_scenarios[dispatch_series[0]])

    #=========================================================================
    ### Second Contour: Difference between Optimum and Suboptimum #2
    
    # find difference between optimal Dispatch type and next available
    diff_2 = performance_metric[Opt_D] - performance_metric[dispatch_series[1]]
    diff_2 /= performance_metric[Opt_D]
    
    # find min and max values of difference
    diff2min = diff_2.min()
    diff2max = diff_2.max()
    
    # define contour levels
    if np.sign(diff2min) != np.sign(diff2max):
        # difference has only one sign, so optimal is strictly better
        levels_2 = [diff2min, 0.7*diff2min, 0.4*diff2min, 0, 0.4*diff2max, 0.9*diff2max, diff2max]
    else:
        # difference changes signs, so sometimes the Dispatch type behaves worse
        levels_2 = [diff2min, 0.75*diff2min, 0.3*diff2min, diff2max]

    # contour plot for Suboptimum #2
    CS2 = axes[n].contour(x_series, y_series, diff_2.T, levels_2, linestyles='--', cmap=cmap_contour)

    # labels for each contour level
    fmt_2 = {}
    for l,s in zip(CS2.levels, levels_2):
        # plotting as percentages
        fmt_2[l] = ' {0:.1f}'.format(s*100) + ' % '
        
    # Add contour labels
    axes[n].clabel(CS2,  inline=True, fmt=fmt_2, fontsize=10)
    
    # creating label for entire contour plot
    CS2.collections[0].set_label('Pct wrt ' + dispatch_scenarios[dispatch_series[1]])
    
    # adding legend for the contours
    axes[n].legend(loc='best')
    
    
