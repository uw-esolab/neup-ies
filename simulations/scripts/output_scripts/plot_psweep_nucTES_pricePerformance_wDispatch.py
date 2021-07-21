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
filename   = 'pricePerfvsDispatch_sizingTESandCycle_irr11pct.nuctes' 
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
npv_aftertax = Storage['npv_aftertax']
irr          = Storage['irr_aftertax']
irr_mean     = np.mean(irr)
    
# =============================================================================
# Set Up Figures and Initialize Lists
# =============================================================================

def print_string_header(header):
    print("\n\n=====================================================")
    print(header)
    print("=====================================================\n")

# creating figure
fig = plt.figure(figsize=[18.5,6])
fig.suptitle("IRR: {0:.1f}".format(irr_mean), fontsize=14, fontweight='bold')
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# adjusting space between 
plt.subplots_adjust(wspace=0.4)

# names of dispatch scenarios
dispatch_scenarios  = ['Pyomo - 24hr Horizon', 'Pyomo - 48hr Horizon', 'No Pyomo - SSC Only']  

# specific performance metrics to plot 
performance_metrics = [annual_energy_array, ppa_array, npv_aftertax]
N              = len(performance_metrics)
get_max        = [True, False, True] # whether each performance is optimized with max or min
axes           = [ax1,ax2,ax3] # list of axes
heat_cmap_list = ['Blues', 'Blues_r', 'Blues_r'] # colormaps used for each performance metric
metric_labels  = ['Annual Energy Production \n(TWh)',
                  'PPA Price \n(cents/kWh)',
                  'NPV \n($M)' ] # labels for each performance metric

# some other arrays
P_array  = P_ref * p_mult #actual power cycle reference outputs
lp = 12 #labelpad

# ========================================   
# Looping!
# ========================================     
print("\n\nBlue -> Better Performance in the 2D Heatmaps")
print("Red  -> Better Performance in the Comparison Contours")
for n in range(N):

    # indeces for dispatch scenarios
    dispatch_series  = list( range(len(dispatch_scenarios) ) ) # declare this every loop
    
    #extracting specific performance metric
    performance_metric = performance_metrics[n]
    
    # get optimum value of performance metric 
    Opt = np.amax(performance_metric) if get_max[n] \
            else np.amin(performance_metric)

    # indeces for this maximum value - { D=Dispatch_Type, T=tshours, P=power_cycle }
    [Opt_D, Opt_T, Opt_P] = [np.where(performance_metric == Opt)[d][0] for d in range(len(performance_metric.shape))]
    
    #=========================================================================
    ### Printing Optimal Values
    optimal_scenario = dispatch_scenarios[Opt_D]
    optimal_metric_value = performance_metric[Opt_D, Opt_T, Opt_P]
    metric_label     = metric_labels[n]
    metric_print_lbl = ''.join(metric_label.split()[:-1])
    metric_name      = metric_label.split()[0]
    metric_unit      = metric_label.split()[-1]
    optimal_tes      = tshours[Opt_T]
    optimal_Pref     = P_array[Opt_P]
    
    prt_opt_header = "Singular Optimal Value happens for: [{0}]".format(optimal_scenario)
    prt_opt_value  = "---- Optimal {0} {1} = {2}".format(metric_print_lbl, metric_unit, optimal_metric_value)
    prt_opt_TES    = "-------- @ tshours = {0}".format(optimal_tes) 
    prt_opt_Pref   = "-------- @ P_ref   = {0}".format(optimal_Pref)
    
    print_string_header("  Plot #{0} - {1}".format(n,metric_print_lbl) )
    print( prt_opt_header )
    print( prt_opt_value  )
    print( prt_opt_TES )
    print( prt_opt_Pref + '\n\n' )
    
    #=========================================================================
    ### Heat Map of Performance Metric
    
    # plotting 2D heat map of the **Dispatch type** that leads to optimum value
    im = axes[n].imshow(performance_metric[Opt_D].T, origin='lower', 
                        cmap=heat_cmap_list[n])

    # setting axis labels
    axes[n].set_xlabel('tshours \n(hr)', fontweight='bold')
    if n == 0:
        axes[n].set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')

    # set title of plot
    title = axes[n].set_title('Optimum: ' + optimal_scenario, fontweight='bold')

    # creating colorbar for the 2D heatmap with label
    cb = fig.colorbar(im, ax=axes[n], fraction=0.06, pad=0.01)
    cb.set_label(metric_label, labelpad= lp, fontweight = 'bold')

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
    ### Contour Maps plotting percent differences
    ### ----- showing what regions of the search space over- and under-perform comparatively
    
    # remove optimal index
    rm = dispatch_series.pop(Opt_D) 
    
    # No Pyomo is typically last in series (if it is not Opt_D)
    comparison_ind      = dispatch_series[0] if Opt_D == 2 else dispatch_series[-1]
    comparison_scenario = dispatch_scenarios[comparison_ind] # str
    comparison_metric   = performance_metric[comparison_ind] # values (2-D)

    # find difference between optimal Dispatch type and next available
    diff     = performance_metric[Opt_D] - comparison_metric
    diff_pct = diff / performance_metric[Opt_D]
    
    # find min and max values of difference pct
    diffmin_pct = diff_pct.min()
    diffmax_pct = diff_pct.max()
    
    # indeces for these max and mins - { T=tshours, P=power_cycle }
    [min_T, min_P] = [np.where(diff_pct == diffmin_pct)[d][0] for d in range(2)]
    [max_T, max_P] = [np.where(diff_pct == diffmax_pct)[d][0] for d in range(2)]
    
    
    # define contour levels
    if np.sign(diffmin_pct) == np.sign(diffmax_pct):
        # difference has only one sign, so optimal is strictly better
        cmap_contour = cm.Reds_r # colormap for contours
        strictly_better = True
        if Opt_D == 2:
            levels = [diffmin_pct, 0.75*diffmin_pct, 0.3*diffmin_pct, diffmax_pct]
        else:
            levels = np.linspace(diffmin_pct, diffmax_pct, 5).tolist()
    else:
        # difference changes signs, so sometimes the Dispatch type behaves worse
        cmap_contour = cm.RdBu_r  if get_max[n] and Opt_D != 2 else cm.RdBu # colormap for contours
        strictly_better = False
        if Opt_D == 2:
            levels = [diffmin_pct, 0.7*diffmin_pct, 0.4*diffmin_pct, 0, 0.4*diffmax_pct, 0.9*diffmax_pct, diffmax_pct]
        else:
            levels = np.hstack([ np.linspace(diffmin_pct, 0, 4) , np.linspace(0, diffmax_pct, 4)[1:] ] )
            levels = levels.tolist()

    # contour plot for Suboptimum #1
    CS = axes[n].contour(x_series, y_series, diff_pct.T, levels, linestyles='--',  cmap=cmap_contour)

    # labels for each contour level
    fmt = {}
    for l,s in zip(CS.levels, levels):
        # plotting as percentages
        fmt[l] = ' {0:.3f}'.format(s*100) + ' % '
        
    # Add contour labels
    axes[n].clabel(CS,  inline=True, fmt=fmt, fontsize=10)
    
    # creating label for entire contour plot
    CS.collections[n].set_label('% Diff w.r.t ' + comparison_scenario)
    
    # adding legend for the contours
    axes[n].legend(loc='best')
        
    #=============================================================================
    ### Printing Outputs

    # are we maximizing or minimizing?
    type_of_opt  = "maximize"   if get_max[n] else "minimize"
    
    # print statements
    prt_comp_scenario   = "Comparison with: {0}".format(comparison_scenario)
    prt_comp_metric     = "(Recall that we want to {0} {1} {2})".format(type_of_opt, metric_name, metric_unit)
    
    print( prt_comp_scenario )
    print( prt_comp_metric + '\n'  )
    
    #----best values
    best_diff      = diff.max() if get_max[n] else diff.min()
    sign_best_diff = "+" if np.sign(best_diff)==1 else "-"
    best_pct       = diff_pct.max() if get_max[n] else diff_pct.min()
    sign_best_pct  = "+" if np.sign(best_pct)==1 else "-"
    tes_best       = tshours[max_T] if get_max[n] else tshours[min_T]
    Pref_best      = P_array[max_P] if get_max[n] else P_array[min_P]

    best_diff = np.abs(best_diff) # already extracted sign out
    best_pct  = np.abs(best_pct)*100
    
    # best print statements
    prt_comp_best       = "++++ [{0}] outperforms [{1}]:".format(optimal_scenario, comparison_scenario)
    prt_comp_best_value = "--------- at best by: {0}{1:.2f} {2} or {3}{4:.4f} % better".format( \
                                                    sign_best_diff, best_diff, metric_unit, sign_best_pct, best_pct) 
    prt_comp_best_tes   = "--------- @ tshours = {0}".format(tes_best) 
    prt_comp_best_Pref  = "--------- @ Pref    = {0}".format(Pref_best) 

    print( prt_comp_best  ) 
    print( prt_comp_best_value )
    print( prt_comp_best_tes )
    print( prt_comp_best_Pref + '\n')
    
    #----worst values
    worst_diff      = diff.min() if get_max[n] else diff.max()
    sign_worst_diff = "+" if np.sign(worst_diff)==1 else "-"
    worst_pct       = diff_pct.min() if get_max[n] else diff_pct.max()
    sign_worst_pct  = "+" if np.sign(worst_pct)==1 else "-"
    tes_worst       = tshours[min_T] if get_max[n] else tshours[max_T]
    Pref_worst      = P_array[min_P] if get_max[n] else P_array[max_P]

    worst_diff = np.abs(worst_diff) # already extracted sign out
    worst_pct  = np.abs(worst_pct)*100
    
    prt_comp_worst       = "++++ [{0}] underperforms relative to  [{1}]:".format(optimal_scenario, comparison_scenario)
    prt_comp_worst_value = "--------- at worst by: {0}{1:.2f} {2} or {3}{4:.4f} % worse".format( \
                                                    sign_worst_diff, worst_diff, metric_unit, sign_worst_pct, worst_pct) 
    prt_comp_worst_tes   = "--------- @ tshours = {0}".format(tes_worst) 
    prt_comp_worst_Pref  = "--------- @ Pref    = {0}".format(Pref_worst) 

    print( prt_comp_worst  ) 
    print( prt_comp_worst_value )
    print( prt_comp_worst_tes )
    print( prt_comp_worst_Pref )