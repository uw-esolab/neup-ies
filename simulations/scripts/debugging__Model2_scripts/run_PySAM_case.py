#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 11:23:19 2022

@author: gabrielsoto
"""

import os, sys
sys.path.append('..')
import modules.DualPlantTES as DualPlantTES
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pint
import numpy as np
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
json = "model2_Hamilton_560_TwinPeaks_x1"  # model2_Hamilton_560_tariffx1 # model2_Hamilton_560_TwinPeaks_x1 # model2_CAISO_Hamilton
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 800  # (MW)
tshours = 2    # (hr)
qdotrec = 500  # (MW)

# defining NE2 module
dptes = DualPlantTES.DualPlantTES(json_name=json, is_dispatch=dispatch, log_dispatch_targets=True) 

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

print("Pref    = {0}".format(dptes.SSC_dict['P_ref'] ) )
print("tshours = {0}".format(dptes.SSC_dict['tshours'] ) )
print("QdotRec = {0}".format(dptes.SSC_dict['q_dot_rec_des'] ) )
# ========================

dptes.run_sim( run_loop=run_loop, export=False, overwrite_dispatch_targets=True )
dt = dptes.Plant
so = dptes.SO

print('Made it past execute.')

# =============================================================================
#     SSC Plots
# =============================================================================

from util.DualPlots import DualPlots
upl = DualPlots(dptes, legend_offset = True, x_shrink=0.7)

# 48 hour plot
fig = plt.figure(figsize=[14,6])

gsQ  = gridspec.GridSpec(3,1, figure=fig)

ax1 = fig.add_subplot( gsQ[0,0], facecolor='white')
ax2 = fig.add_subplot( gsQ[1:,0], facecolor='white')

plt_allTime = False
title = 'SSC Results - 48 hrs'
start = 24 * 0
end   = 24 * 364

upl.plot_SSC_power_and_energy(ax1 , plot_all_time=plt_allTime, title_label=title, start_hr=start, end_hr=end, hide_x=True, x_legend=1.2, y_legend_L=1.0, y_legend_R=0.3)
upl.plot_SSC_op_modes(ax2, plot_all_time=plt_allTime, start_hr=start, end_hr=end, hide_x=False )
# upl.plot_SSC_massflow(ax3, plot_all_time=plt_allTime, start_hr=start, end_hr=end, y_legend_L=0.8, y_legend_R=0.3)


# =============================================================================
#   Creating Pyomo Plotting Object
# =============================================================================

# retrieving the DispatchPlots class
from util.DualPlots import DualDispatchPlots
# specifying dispatch model
ind = 11
# extracting specific, solved dispatch model
dm = dptes.disp_models[str(ind)]
# create Dispatch plotting object
dpl = DualDispatchPlots(dm, lp=8, legend_offset=True, x_shrink=0.75)


# =============================================================================
#   Pyomo Plotting
# =============================================================================

fig = plt.figure(figsize=[12, 10])
ax1 = fig.add_subplot(511)
ax2 = fig.add_subplot(512)
# ax3 = fig.add_subplot(513)
ax4 = fig.add_subplot(514)
ax5 = fig.add_subplot(515)

plt.subplots_adjust(hspace=0)
# plt.gcf().subplots_adjust(bottom=0.1) # leaving room at the bottom

dpl.plot_pyomo_energy(      ax1, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_power(       ax2, x_legend=1.15, y_legend_L=1.0, hide_x=True )
# dpl.plot_pyomo_power_ramps( ax3, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_dual_bin(    ax4, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_cycle_bin(   ax5, x_legend=1.04, y_legend_L=1.0)

fig.suptitle('Pyomo Optimization Start at Day {0}'.format(ind), fontsize=16)