#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 12:08:24 2022

@author: gabrielsoto
"""

import os,sys,copy
sys.path.append('..')
import modules.DualIndirectTES as DualIndirectTES
import matplotlib.pyplot as plt
import pint
import numpy as np
import pyomo.environ as pe
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
#     Run Simulation
# =============================================================================

# modifying inputs
json = "model2b_Hamilton_FS_tariffx1"  # model2b_Hamilton_FS_TwinPeaks_x1 # model2b_Hamilton_FS_tariffx1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 700  # (MW)
tshours = 8    # (hr)
qdotrec = 500  # (MW)

# defining NE2 module
dptes = DualIndirectTES.DualIndirectTES(json_name=json, is_dispatch=dispatch, log_dispatch_targets=True) 

# horizons
dptes.PySAM_dict['ssc_horizon']   = sscH
dptes.ssc_horizon   = sscH * dptes.u.hr
dptes.PySAM_dict['pyomo_horizon'] = pyoH
dptes.pyomo_horizon = pyoH * dptes.u.hr

# saving/updating PYSAM dict to nuctes
dptes.dispatch_wrap = dptes.create_dispatch_wrapper( dptes.PySAM_dict )

# update SSC dictionary parameters
dptes.SSC_dict['P_ref'] = Pref
dptes.SSC_dict['tshours'] = tshours
dptes.SSC_dict['q_dot_rec_des'] = qdotrec

dptes.dispatch_wrap.set_design()

print(json )
print("Pref    = {0}".format(dptes.SSC_dict['P_ref'] ) )
print("tshours = {0}".format(dptes.SSC_dict['tshours'] ) )
print("QdotRec = {0}".format(dptes.SSC_dict['q_dot_rec_des'] ) )
# ========================

dptes.run_sim( run_loop=run_loop, export=False, overwrite_dispatch_targets=True )
dt = dptes.Plant
so = dptes.SO

print('Made it past execute.')

# =============================================================================
#   Creating Pyomo Plotting Object
# =============================================================================

# specifying dispatch model
ind = 0
# extracting specific, solved dispatch model
dm = dptes.disp_models[str(ind)]

# lambda functions
extract_from_model = lambda name: getattr(dm.model, name)
extract_array      = lambda name: np.array([ pe.value(extract_from_model(name)[t]) for t in dm.model.T ])
extract_energy     = lambda name: (extract_array(name)*u.kWh).to('MWh') 
extract_power      = lambda name: (extract_array(name)*u.kW).to('MW') 

# Time Array
t_full = extract_array('Delta_e') * u.hr

# Power Arrays
wdot_array    = extract_power('wdot').m
xnp_array     = extract_power('xnp').m
xntes_array   = extract_power('xntes').m
xtesp_array   = extract_power('xtesp').m
xr_array      = extract_power('xr').m

# Energy In Tank Array
s_array   = extract_energy('s').m

# Binary Arrays
ytesp_array   = extract_array('ytesp')
# yntes_array   = extract_array('yntes')


# =============================================================================
# plotting efficiencies
# =============================================================================

# Power Cycle Parameters
Wnc = pe.value( extract_from_model('Wnc') ) / 1000
Wu = pe.value( extract_from_model('Wdotu') )
Qu = pe.value( extract_from_model('Qu') )
etaLD = pe.value( extract_from_model('eta_LD') )
etaHD = pe.value( extract_from_model('eta_HD') )
eta = pe.value( extract_from_model('eta_des') )

b1 = Wu - etaLD*Qu
b2 = Wu - etaHD*Qu

x = np.linspace(0, Pref / 0.51, 1000) * 1000
yLD = etaLD*x + b1
yHD = etaHD*x + b2

x /= 1000
yLD /= 1000
yHD /= 1000




# Thermal Power Array
qdot_array = copy.deepcopy(wdot_array)
# qdot_array /= etaLD
qdot_array[wdot_array < Wnc] /= eta
qdot_array[wdot_array > Wnc] /= eta

sorted_q_ind = np.argsort(xnp_array + xtesp_array) 
q_dot_sorted = xnp_array + xtesp_array
q_dot_sorted = q_dot_sorted[sorted_q_ind]
w_dot_sorted = wdot_array[sorted_q_ind]

# fig = plt.figure()
# ax = fig.gca()
# ax.plot( x, yLD, 'C1--', label='Low Demand Efficiency')
# ax.plot( x, yHD, 'C2--', label='High Demand Efficiency')
# ax.plot(q_dot_sorted, w_dot_sorted, 'o')

        
# =============================================================================
#   Pyomo Plotting
# =============================================================================
t_full = t_full.m

def plot_with_marker(ax, parray, color, linewidth, label ):
    ax.plot(t_full, parray, color=color, linewidth=linewidth, label=label)
    ax.plot(t_full, parray, color=color, marker='o')
    

lw = 3
fig = plt.figure(figsize=[12, 10])
ax = fig.add_subplot(311)

ax.set_title( f'Pref = {Pref} MWe, \n Tshours = {tshours}, \n CSP Output = {qdotrec}, \n json = {json}' )

plot_with_marker(ax, qdot_array, color='C3', linewidth=lw, label='PC Power Out (MWt)') 
plot_with_marker(ax, xnp_array+xtesp_array, color='k', linewidth=lw,    label='Total Est. Power In (MWt)') 
plot_with_marker(ax, xnp_array, color='darkblue',  linewidth=lw, label='LFR Power In (MWt)')
plot_with_marker(ax, xtesp_array, color='teal', linewidth=lw,    label='TES Power In (MWt)') 
ax.set_ylabel('Power Cycle CV', fontweight='bold')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
ax.grid(True)

ax = fig.add_subplot(312)
plot_with_marker(ax, xnp_array+xntes_array, color='black',  linewidth=lw, label='Total Power (MWt)') 
plot_with_marker(ax, xnp_array, color='darkblue',  linewidth=lw,          label='LFR Power to PC (MWt)')
plot_with_marker(ax, xntes_array, color='darkorange',  linewidth=lw,      label='LFR Power to TES (MWt)') 
plot_with_marker(ax, xr_array, color='darkred',  linewidth=lw,      label='CSP Power to TES (MWt)') 
ax.set_ylabel('Nuclear CV', fontweight='bold')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
ax.grid(True)

ax = fig.add_subplot(313)
plot_with_marker(ax, s_array, color='purple',  linewidth=lw, label='') 
ax.set_ylabel('TES Charge State', fontweight='bold')
ax.set_xlabel('Time (hrs)', fontweight='bold')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

