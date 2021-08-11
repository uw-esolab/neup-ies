#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 13:17:59 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
import os, pint, time, copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# set up
# =============================================================================

# array of successful P_cycle and tshours combinations
# p_cycle_0     = np.array([ 700, 650, 600, 550,      450, 400, 350, 300, 250, 200, 150 ]) # Pref= 500, tshours=0
# p_cycle_2     = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=2
# p_cycle_4     = np.array([ 700, 650, 600, 550, 500,                300, 250, 200, 150 ]) # Pref= [350,450], tshours=4
# p_cycle_6     = np.array([ 700, 650, 600, 550, 500,      400,      300, 250, 200, 150 ]) # Pref= {350,450}, tshours=6
# p_cycle_8     = np.array([ 700, 650, 600, 550, 500,           350, 300, 250,      150 ]) # Pref= {200,400,450}, tshours=8
# p_cycle_10    = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=10
# p_cycle_12    = np.array([ 700, 650, 600, 550, 500,           350, 300, 250,      150 ]) # Pref= {200,400,450}, tshours=12
# p_cycle_14    = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=14

# setting up arrays to cycle through
tshours    = np.array([ 0, 2, 4, 6, 8, 10, 12, 14 ])
p_cycle    = np.array([ 1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150 ]) 

# exceptions!
exceptions = {
    "0":  [500],
    "2":  [450],
    "4":  [350, 400, 450],
    "6":  [350, 450],
    "8":  [200, 400, 450],
    "10": [450],
    "12": [200, 400, 450], 
    "14": [450]
    }

# formally defining the iterator arrays
iterator1 = tshours  
iterator2 = p_cycle

# initializing an empty array
empty = np.zeros([  len(iterator1) , len(iterator2) ])
empty2D = np.zeros([  len(iterator1) , len(iterator2), 2 ])

# initializing output arrays
defocus  = empty2D.copy()
gen      = empty2D.copy()
qdot_nuc = empty2D.copy()
TES_CH = empty.copy()
TES_DC = empty.copy()

# =============================================================================
# Operating Modes analysis
# =============================================================================

from util.PostProcessing import Plots
tmp_modes = lambda: None
Plots.set_operating_modes_list(tmp_modes)

op_modes_list = tmp_modes.operating_modes


# =============================================================================
# Double Loop
# =============================================================================

# starting the time counter
tic = time.perf_counter()

# TRIPLE LOOP
for i, th in enumerate(iterator1): #over tshours
    for j, pc in enumerate(iterator2): #over P_cycle rating
        
        idx = (i,j)
        mean_idx = (i,j,0)
        stdv_idx = (i,j,1)

        # print current position in loop
        print("output tshours :   ", th)
        print("output Pref    :   ", pc)
        
        # initialize the PySSC Wrapper
        pw = PySSCWrapper(json_name='model1',is_debug=True)
        
        # update SSC dictionary parameters
        pw.sscdict['tshours'] = float(th)
        pw.sscdict['P_ref']   = float(pc)
        
        if pc in exceptions[str(th)]:
            dummy_val = -1
            # log output means
            gen[mean_idx]       = dummy_val
            qdot_nuc[mean_idx]  = dummy_val
            defocus[mean_idx]   = dummy_val
            
            # log output std dev
            gen[stdv_idx]       = dummy_val
            qdot_nuc[stdv_idx]  = dummy_val
            defocus[stdv_idx]   = dummy_val
            
            TES_CH[idx] = dummy_val
            TES_DC[idx] = dummy_val
            
        else:
            # run simulationstd_label
            success = True
            try:
                # run SSC through PySSC
                pw.run_sim()
            except Exception as err:
                print(" Run failed.")
                success = False
            
            if success:
                
                # analyzing OP modes
                op_mode_array       = pw.get_array('op_mode_1') 
                op_mode_str_profile = [op_modes_list[ int(x) ] for x in op_mode_array]
                
                TES_charging    = np.array(["TES_CH" in x for x in op_mode_str_profile])
                TES_discharging = np.array(["TES_DC" in x for x in op_mode_str_profile])
                
                TES_CH[idx] = int(TES_charging.sum() > 0)
                TES_DC[idx] = int(TES_discharging.sum() > 0)
                
                # log output means
                gen[mean_idx]       = np.mean( pw.get_array('gen') / 1e3     )
                qdot_nuc[mean_idx]  = np.mean( pw.get_array('q_dot_rec_inc') )
                defocus[mean_idx]   = np.mean( pw.get_array('defocus') )
                
                # log output std dev
                gen[stdv_idx]       = np.std( pw.get_array('gen') / 1e3     )
                qdot_nuc[stdv_idx]  = np.std( pw.get_array('q_dot_rec_inc') )
                defocus[stdv_idx]   = np.std( pw.get_array('defocus') )
            
        # reset the Plant and Grid, prevents memory leak
        del pw
            
# end time counter
toc = time.perf_counter()        

# =============================================================================
# Plot
# =============================================================================
cmap = cm.hot_r

array = copy.deepcopy( qdot_nuc )
title_label = "LFR_Qdot" #"'Defocus' Fraction"
mean_label = "LFR_Qdot Mean \n(MW)"
std_label = "LFR_Qdot Std Dev \n(MW)"
# ========== Arrays ==========

# fixing "bad" indeces
min_array_avg = np.abs(array[:,:,0]).min()
fix0,fix1       = np.where( array[:,:,0] == -1)
defocus[fix0,fix1,0] = -min_array_avg

# fixing "bad" indeces
min_array_std = np.abs(array[:,:,1]).min()
fix0,fix1       = np.where( array[:,:,1] == -1)
array[fix0,fix1,1] = -min_array_std


# create figure
fig = plt.figure()
ax1  = fig.add_subplot(121)
ax2  = fig.add_subplot(122)
fig.suptitle(title_label, fontweight='bold')

asp_df = 0.7
im1 = ax1.imshow(array[:,:,0].T, origin='upper', cmap=cmap, aspect=asp_df)
im2 = ax2.imshow(array[:,:,1].T, origin='upper', cmap=cmap, aspect=asp_df)

# setting axis labels
ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
ax2.set_xlabel('tshours\n(hr)', fontweight='bold')

ax1.set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')


# creating colorbar for the 2D heatmap with label
cb1 = fig.colorbar(im1, ax=ax1, fraction=0.06, pad=0.01)
cb1.set_label(mean_label, labelpad= 8, fontweight = 'bold')

# creating colorbar for the 2D heatmap with label
cb2 = fig.colorbar(im2, ax=ax2, fraction=0.06, pad=0.01)
cb2.set_label(std_label, labelpad= 8, fontweight = 'bold')

# setting tick marks for x and y axes
ax1.set_xticks(range(len(tshours)))
ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )

ax2.set_xticks(range(len(tshours)))
ax2.set_xticklabels( ['{0}'.format(t) for t in tshours] )

ax1.set_yticks(range(len(p_cycle)))
ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle] )
ax2.axes.yaxis.set_visible(False)

# =============================================================================
# Boolean heatmaps
# =============================================================================
cmap = cm.hot_r

array1 = copy.deepcopy( TES_CH )
array2 = copy.deepcopy( TES_DC )

array1_label = "Did TES Charge?"
array2_label = "Did TES Discharge?"
# ========== Arrays ==========


# create figure
fig = plt.figure()
ax1  = fig.add_subplot(121)
ax2  = fig.add_subplot(122)
fig.suptitle("PySSC run - No Pyomo", fontweight='bold')

asp_df = 0.7
im1 = ax1.imshow(array1.T, origin='upper', cmap=cmap, aspect=asp_df)
im2 = ax2.imshow(array2.T, origin='upper', cmap=cmap, aspect=asp_df)

# setting axis labels
ax1.set_xlabel('tshours\n(hr)', fontweight='bold')
ax2.set_xlabel('tshours\n(hr)', fontweight='bold')

ax1.set_ylabel('Power Cycle Output\n(MW)', fontweight='bold')


# creating colorbar for the 2D heatmap with label
cb1 = fig.colorbar(im1, ax=ax1, fraction=0.06, pad=0.01)
cb1.set_label(array1_label, labelpad= 8, fontweight = 'bold')

cb1.set_ticks([-1,0,1])
cb1.set_ticklabels( ["N/A", "No", "Yes"] )

# creating colorbar for the 2D heatmap with label
cb2 = fig.colorbar(im2, ax=ax2, fraction=0.06, pad=0.01)
cb2.set_label(array2_label, labelpad= 8, fontweight = 'bold')

cb2.set_ticks([-1,0,1])
cb2.set_ticklabels( ["N/A", "No", "Yes"] )

# setting tick marks for x and y axes
ax1.set_xticks(range(len(tshours)))
ax1.set_xticklabels( ['{0}'.format(t) for t in tshours] )

ax2.set_xticks(range(len(tshours)))
ax2.set_xticklabels( ['{0}'.format(t) for t in tshours] )

ax1.set_yticks(range(len(p_cycle)))
ax1.set_yticklabels( ['{0:.0f}'.format(P) for P in p_cycle] )
ax2.axes.yaxis.set_visible(False)