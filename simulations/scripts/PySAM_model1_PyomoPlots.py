#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:31:26 2021

@author: gabrielsoto
"""

import modules.NuclearTES as NuclearTES
import os
import matplotlib.pyplot as plt
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
#   Run Simulation
# =============================================================================

# creating NE2 module object
nuctes = NuclearTES.NuclearTES(is_dispatch=True)
# run simulation
nuctes.run_sim( run_loop=True )
# extract solved module
nt = nuctes.Plant


# =============================================================================
#   Creating Pyomo Plotting Object
# =============================================================================

# retrieving the DispatchPlots class
from util.PostProcessing import DispatchPlots
# specifying dispatch model
ind = 100
# extracting specific, solved dispatch model
dm = nuctes.disp_models[str(ind)]
# create Dispatch plotting object
dpl = DispatchPlots(dm, lp=8, legend_offset=True, x_shrink=0.75)


# =============================================================================
#   Pyomo Plotting
# =============================================================================

fig = plt.figure(figsize=[12, 10])
ax1 = fig.add_subplot(511)
ax2 = fig.add_subplot(512)
ax3 = fig.add_subplot(513)
ax4 = fig.add_subplot(514)
ax5 = fig.add_subplot(515)

plt.subplots_adjust(hspace=0)
# plt.gcf().subplots_adjust(bottom=0.1) # leaving room at the bottom

dpl.plot_pyomo_energy(      ax1, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_power(       ax2, x_legend=1.15, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_power_ramps( ax3, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_nuclear_bin( ax4, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_cycle_bin(   ax5, x_legend=1.04, y_legend_L=1.0)

