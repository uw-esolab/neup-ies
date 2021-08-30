#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 11:27:41 2021

@author: gabrielsoto
"""


from util.PySSCWrapper import PySSCWrapper
import os, pint, time
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

# setting up arrays to cycle through
# p_cycle_0     = np.array([ 700, 650, 600, 550,      450, 400, 350, 300, 250, 200, 150 ]) # Pref= 500, tshours=0
# p_cycle_2     = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=2
# p_cycle_4     = np.array([ 700, 650, 600, 550, 500,                300, 250, 200, 150 ]) # Pref= [350,450], tshours=4
# p_cycle_6     = np.array([ 700, 650, 600, 550, 500,      400,      300, 250, 200, 150 ]) # Pref= {350,450}, tshours=6
p_cycle_8     = np.array([ 700, 650, 600, 550, 500,           350, 300, 250,      150 ]) # Pref= {250,400,450}, tshours=8
# p_cycle_10    = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=10
# p_cycle_12    = np.array([ 700, 650, 600, 550, 500,           350, 300, 250,      150 ]) # Pref= {200,400,450}, tshours=12
# p_cycle_14    = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=14
# p_cycle    = np.array([ 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150 ]) 
tshours = 8

# formally defining the iterator arrays
iterator1 = p_cycle_8 

# initializing an empty array
year_in_hours = int(  8760  )
empty = np.zeros([  len(iterator1) , year_in_hours ])

# initializing output arrays
time_array     = empty.copy()
defocus_array  = empty.copy()
gen_array      = empty.copy()
eta_array      = empty.copy()
qdot_nuc_array = empty.copy()
op_mode_array  = empty.copy()
sim_success    = empty.copy()


# starting the time counter
tic = time.perf_counter()

# TRIPLE LOOP
for i, pc in enumerate(iterator1): #over P_cycle rating

    # print current position in loop
    print("output Pref :   ", pc)
    
    # initialize the PySSC Wrapper
    pw = PySSCWrapper(json_name='model1',is_debug=True)
    
    # update SSC dictionary parameters
    pw.sscdict['P_ref']   = float(pc)
    pw.sscdict['tshours'] = float(tshours)
    
    # run simulation
    success = True
    try:
        # run SSC through PySSC
        pw.run_sim()
    except Exception as err:
        print(" Run failed.")
        success = False
    
    if success:
        # log outputs
        time_array[i,:]      = pw.get_array('time_hr') / 24
        gen_array[i,:]       = pw.get_array('gen') / 1e3
        eta_array[i,:]       = pw.get_array('eta')
        qdot_nuc_array[i,:]  = pw.get_array('q_dot_rec_inc')
        defocus_array[i,:]   = pw.get_array('defocus')
        op_mode_array[i,:]   = pw.get_array('op_mode_1')

    sim_success[i,:] = True
    
    # reset the Plant and Grid, prevents memory leak
    del pw
            
# end time counter
toc = time.perf_counter()        


# =============================================================================
# Plot
# =============================================================================
cmap = cm.viridis if iterator1[-1] > iterator1[0] else cm.viridis_r

fig = plt.figure(figsize=(10,6))
ax  = fig.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

N = len(iterator1)
for j, p in enumerate(iterator1):
    j_value = j / N 
    c_value = cmap(j_value)
    
    ax.plot(time_array[j,:], defocus_array[j,:], linewidth=2, \
            color=c_value, label="P_ref = {0} MW".format(p))
    
ax.set_xlabel("Time (hr)", fontweight="bold")
ax.set_ylabel("Defocus", fontweight="bold")
ax.set_title("Default tariff schedule \ntshours = {0} hr".format(tshours), fontweight="bold")
ax.legend(  bbox_to_anchor=(1, 0.8)) 

# eta plot
fig = plt.figure(figsize=(10,6))
ax  = fig.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

N = len(iterator1)
for j, p in enumerate(iterator1):
    j_value = j / N 
    c_value = cmap(j_value)
    
    ax.plot(time_array[j,:], eta_array[j,:], linewidth=2, \
            color=c_value, label="P_ref = {0} MW".format(p))
    
ax.set_xlabel("Time (hr)", fontweight="bold")
ax.set_ylabel("Eta", fontweight="bold")
ax.set_title("Eta \ntshours = {0} hr".format(tshours), fontweight="bold")
ax.legend(  bbox_to_anchor=(1, 0.8)) 


# power plot
fig = plt.figure(figsize=(10,6))
ax  = fig.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

N = len(iterator1)
for j, p in enumerate(iterator1):
    j_value = j / N 
    c_value = cmap(j_value)
    
    ax.plot(time_array[j,:], gen_array[j,:], linewidth=2, \
            color=c_value, label="P_ref = {0} MW".format(p))
    
ax.set_xlabel("Time (hr)", fontweight="bold")
ax.set_ylabel("Power Generated (MW)", fontweight="bold")
ax.set_title("Power Generated \ntshours = {0} hr".format(tshours), fontweight="bold")
ax.legend(  bbox_to_anchor=(1, 0.8)) 

