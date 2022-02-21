#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 12:20:12 2021

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import modules.NuclearTES as NuclearTES
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
json = "model1_CAISO"   # model1_CAISO # model1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 700  # (MW)
tshours = 10    # (hr)

# ========================

# defining directories
nuctes = NuclearTES.NuclearTES(json_name=json, is_dispatch=dispatch, log_dispatch_targets=False) 
output_file = 'output.csv'

# horizons
nuctes.PySAM_dict['ssc_horizon']   = sscH
nuctes.ssc_horizon   = sscH * nuctes.u.hr
nuctes.PySAM_dict['pyomo_horizon'] = pyoH
nuctes.pyomo_horizon = pyoH * nuctes.u.hr

# saving/updating PYSAM dict to nuctes
nuctes.dispatch_wrap = nuctes.create_dispatch_wrapper( nuctes.PySAM_dict )

# update SSC dictionary parameters
nuctes.SSC_dict['P_ref'] = Pref
nuctes.SSC_dict['tshours'] = tshours
nuctes.dispatch_wrap.set_design()

print("Pref    = {0}".format(nuctes.SSC_dict['P_ref'] ) )
print("tshours = {0}".format(nuctes.SSC_dict['tshours'] ) )
# ========================

nuctes.run_sim( run_loop=run_loop, export=False, filename=output_file )
nt = nuctes.Plant
so = nuctes.SO

print('Made it past execute.')

# =============================================================================
#   Creating Pyomo Plotting Object
# =============================================================================

# retrieving the DispatchPlots class
from util.PostProcessing import DispatchPlots
# specifying dispatch model
ind = 138
# extracting specific, solved dispatch model
dm = nuctes.disp_models[str(ind)]
# create Dispatch plotting object
dpl = DispatchPlots(dm, lp=8, legend_offset=True, x_shrink=0.75)


# =============================================================================
#   Pyomo Plotting
# =============================================================================

fig = plt.figure(figsize=[12, 10])
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)

plt.subplots_adjust(hspace=0)
# plt.gcf().subplots_adjust(bottom=0.1) # leaving room at the bottom

title_pyo = 'Pyomo Results \nTime after Sim Start = {0} d, Pyomo Horizon = {1} hr \nPref = {2:.2f} , tshours = {3}'.format( \
                                (ind*sscH*u.hr).to('d').m, pyoH, nuctes.SSC_dict['P_ref'], nuctes.SSC_dict['tshours'])

dpl.plot_pyomo_energy(      ax1, x_legend=1.04, y_legend_L=1.0, hide_x=True, title_label=title_pyo )
dpl.plot_pyomo_power(       ax2, x_legend=1.15, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_nuclear_bin( ax3, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_cycle_bin(   ax4, x_legend=1.04, y_legend_L=1.0)