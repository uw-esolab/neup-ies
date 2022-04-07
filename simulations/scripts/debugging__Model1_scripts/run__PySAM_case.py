#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 14:53:16 2021

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
json = "model1_CAISO_Hamilton"   # model1_CAISO_Hamilton # model1_Hamilton_560_tariffx2
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 1000 # (MW)
tshours = 1    # (hr)

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

print(json )
print("Pref    = {0}".format(nuctes.SSC_dict['P_ref'] ) )
print("tshours = {0}".format(nuctes.SSC_dict['tshours'] ) )
# ========================

nuctes.run_sim( run_loop=run_loop, export=False, filename=output_file )
nt = nuctes.Plant
so = nuctes.SO

print('Made it past execute.')

# =============================================================================
#     Display Results
# =============================================================================

annual_energy          = (nt.Outputs.annual_energy*u.kWh).to('TWh')
capacity_factor        = nuctes.capacity_factor.magnitude * 100
annual_total_water_use = nt.Outputs.annual_total_water_use
ppa                    = so.Outputs.ppa
lppa_nom               = so.Outputs.lppa_nom
lppa_real              = so.Outputs.lppa_real
lcoe_nom               = so.Outputs.lcoe_nom
lcoe_real              = so.Outputs.lcoe_real
project_return_aftertax_npv = so.Outputs.project_return_aftertax_npv/1e6
flip_actual_irr        = so.Outputs.flip_actual_irr
flip_actual_year       = so.Outputs.flip_actual_year
project_return_aftertax_irr = so.Outputs.project_return_aftertax_irr
cost_installed         = so.Outputs.cost_installed/1e6
size_of_equity         = so.Outputs.size_of_equity/1e6
size_of_debt           = so.Outputs.size_of_debt/1e6

print('')
print('                        Nuclear')
print ('Annual energy (year 1)        =   ', annual_energy)
print ('Capacity factor (year 1)      =   ', capacity_factor, ' %')
print ('Annual Water Usage            =   ', annual_total_water_use, ' m3')
print ('PPA price (year 1)            =   ', ppa, ' ¢/kWh')
print ('Levelized PPA price (nominal) =   ', lppa_nom, ' ¢/kWh')
print ('Levelized PPA price (real)    =   ', lppa_real, ' ¢/kWh')
print ('Levelized COE (nominal)       =   ', lcoe_nom, ' ¢/kWh')
print ('Levelized COE (real)          =   ', lcoe_real, ' ¢/kWh')
print ('Net present value             =  $', project_return_aftertax_npv, ' M')
print ('Internal rate of return (IRR) =   ', flip_actual_irr, ' %')
print ('Year IRR is achieved          =   ', flip_actual_year)
print ('IRR at end of project         =   ', project_return_aftertax_irr, ' %')
print ('Net capital cost              =  $', cost_installed, ' M')
print ('Equity                        =  $', size_of_equity, ' M')
print ('Size of debt                  =  $', size_of_debt, ' M')

# =============================================================================
#     Plotting
# =============================================================================

from util.PostProcessing import Plots
upl = Plots(nuctes, legend_offset = True, x_shrink=0.7, fE_min=0, fE_max=1.05)

# 48 hour plot
fig = plt.figure(figsize=[14,6])
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

plt_allTime = False
title = 'SSC Results - {0}\nSSC horizon = {1} hr, Pyomo Horizon = {2} hr \nPref = {3:.2f} , tshours = {4}'.format( \
                                json, sscH, pyoH, nuctes.SSC_dict['P_ref'], nuctes.SSC_dict['tshours'])
start = 24*0
end   = 24*364

upl.plot_SSC_power_and_energy(ax1 , plot_all_time=plt_allTime, title_label=title, start_hr=start, end_hr=end, hide_x=True, x_legend=1.2, y_legend_L=1.0, y_legend_R=0.3)
upl.plot_SSC_op_modes(ax2, plot_all_time=plt_allTime, start_hr=start, end_hr=end, hide_x=True )
upl.plot_SSC_massflow(ax3, plot_all_time=plt_allTime, start_hr=start, end_hr=end, days_on_x=True, y_legend_L=0.8, y_legend_R=0.3)


# =============================================================================
#   Creating Pyomo Plotting Object
# =============================================================================

# retrieving the DispatchPlots class
from util.PostProcessing import DispatchPlots
# specifying dispatch model
ind = 292
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

title_pyo = 'Pyomo Results \nTime after Sim Start = {0} d, Pyomo Horizon = {1} hr \nPref = {2:.2f} , tshours = {3}'.format( \
                                (ind*sscH*u.hr).to('d').m, pyoH, nuctes.SSC_dict['P_ref'], nuctes.SSC_dict['tshours'])

dpl.plot_pyomo_energy(      ax1, x_legend=1.04, y_legend_L=1.0, hide_x=True, title_label=title_pyo )
dpl.plot_pyomo_power(       ax2, x_legend=1.15, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_power_ramps( ax3, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_nuclear_bin( ax4, x_legend=1.04, y_legend_L=1.0, hide_x=True )
dpl.plot_pyomo_cycle_bin(   ax5, x_legend=1.04, y_legend_L=1.0)
