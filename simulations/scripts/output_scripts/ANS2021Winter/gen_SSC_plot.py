#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 09:15:26 2021

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
json = "model1"   # model1_CAISO # model1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 525  # (MW)
tshours = 8    # (hr)

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
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

plt_allTime = True
title = 'SSC Results \nSSC horizon = {0} hr, Pyomo Horizon = {1} hr \nPref = {2:.2f} , tshours = {3}'.format( \
                                sscH, pyoH, nuctes.SSC_dict['P_ref'], nuctes.SSC_dict['tshours'])
start = 0
end   = start + 3600

upl.plot_SSC_power_and_energy(ax1 , plot_all_time=plt_allTime, title_label=title, start_hr=start, end_hr=end, hide_x=True, x_legend=1.2, y_legend_L=1.0, y_legend_R=0.3)
upl.plot_SSC_op_modes(ax2, plot_all_time=plt_allTime, start_hr=start, end_hr=end)

