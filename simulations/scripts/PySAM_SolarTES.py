#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:13:14 2021

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import modules.SolarTES as SolarTES
import matplotlib.pyplot as plt
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
tshours = 4   # (hr)


# defining NE2 module
soltes = SolarTES.SolarTES( is_dispatch=True )

# update SSC dictionary parameters
# soltes.SSC_dict['P_ref'] = Pref
soltes.SSC_dict['tshours'] = tshours
# soltes.SSC_dict['q_dot_rec_des'] = qdotrec

soltes.dispatch_wrap.set_design()


soltes.run_sim( run_loop=True, export=False )

st = soltes.Plant
so = soltes.SO

print('Made it past execute.')

# =============================================================================
#     Display Results
# =============================================================================

annual_energy          = (st.Outputs.annual_energy*u.kWh).to('TWh')
capacity_factor        = soltes.capacity_factor.magnitude * 100
annual_total_water_use = st.Outputs.annual_total_water_use
ppa                    = so.Outputs.ppa
lppa_nom               = so.Outputs.lppa_nom
lppa_real              = so.Outputs.lppa_real
lcoe_nom               = so.Outputs.lcoe_nom
lcoe_real              = so.Outputs.lcoe_real
project_return_aftertax_npv = so.Outputs.project_return_aftertax_npv / 1e6
flip_actual_irr        = so.Outputs.flip_actual_irr
flip_actual_year       = so.Outputs.flip_actual_year
project_return_aftertax_irr = so.Outputs.project_return_aftertax_irr
cost_installed         = so.Outputs.cost_installed / 1e6
size_of_equity         = so.Outputs.size_of_equity / 1e6
size_of_debt           = so.Outputs.size_of_debt / 1e6

print('')
print('                        CSP')
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

from util.SolarPlots import SolarPlots
upl = SolarPlots(soltes, legend_offset = True, x_shrink=0.7)

# 48 hour plot
fig = plt.figure(figsize=[14,6])
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

plt_allTime = False
title = 'SSC Results - 48 hrs'
start = 24 * 0
end   = 24 * 364

upl.plot_SSC_power_and_energy(ax1 , plot_all_time=plt_allTime, title_label=title, start_hr=start, end_hr=end, hide_x=True, x_legend=1.2, y_legend_L=1.0, y_legend_R=0.3)
upl.plot_SSC_op_modes(ax2, plot_all_time=plt_allTime, start_hr=start, end_hr=end, hide_x=True )
upl.plot_SSC_massflow(ax3, plot_all_time=plt_allTime, start_hr=start, end_hr=end, y_legend_L=0.8, y_legend_R=0.3)