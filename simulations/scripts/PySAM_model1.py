#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:47:58 2020

Most recently tested against PySAM 2.1.4

@author: frohro
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

# defining directories
nuctes = NuclearTES.NuclearTES()
output_file = 'output.csv'
nuctes.run_sim( run_loop=False, export=True, filename=output_file )
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
upl = Plots(nuctes)

# 48 hour plot
fig = plt.figure(figsize=[10,8])
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

plot_full_time = True

upl.plot_SSC_power_and_energy(ax1 )#, plot_full_time=plot_full_time, title_label='SSC Results - 48 hrs')
upl.plot_SSC_op_modes(ax2) #, plot_full_time=plot_full_time)
upl.plot_SSC_massflow(ax3) #, plot_full_time=plot_full_time)

# full 1 year plot
# figF = plt.figure(figsize=[10,8])
# ax1F = figF.add_subplot(311)
# ax2F = figF.add_subplot(312)
# ax3F = figF.add_subplot(313)

# upl.plot_SSC_power_and_energy(ax1F, True, title_label='SSC Results - Full Year')
# upl.plot_SSC_op_modes(ax2F, True)
# upl.plot_SSC_massflow(ax3F, True)


# upl.plot_SSC_op_modes(False)
# upl.plot_SSC_op_modes(True)
# upl.plot_SSC_power_and_energy(False)
# upl.plot_SSC_power_and_energy(True)
# upl.plot_SSC_massflow(False)
# upl.plot_SSC_massflow(True)


# plt.figure()
# plt.plot(t_plot,nuctes.df_array)


# fig = plt.figure(figsize=[10,8])
# ax1 = fig.gca()
# # ax1 = fig.add_subplot(311)
# # ax2 = fig.add_subplot(312)
# # ax3 = fig.add_subplot(313)

# # Energy plot
# ax1.plot(t_plot, e_ch_tes, linewidth = lw, label='Salt Charge Level (MWht)')
# ax1.plot(t_plot, p_cycle, linewidth = lw, label='P_cycle (Electric)')
# ax1.plot(t_plot, q_dot_rec_in, linewidth = lw, label='Q_dot to Salt (Thermal)')
# ax1.plot(t_plot, gen, linewidth = lw, label='Power generated')
# ax1.set_ylabel('Power (MW)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax1.legend(loc=loc,fontsize=fsl)
# ax1.set_xlabel('Time (days)', labelpad=lp, fontsize=fs, fontweight='bold')

# Mass flow rate plot
# ax2.plot(t_plot, m_dot, linewidth = lw, label='m_dot_water_pc')
# ax2.set_ylabel('Mass flow (kg/s)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax2.legend(loc=loc,fontsize=fsl)

# Temperature plot
# ax3.plot(t_plot, T_pc_in, linewidth = lw, label='PC HTF inlet')
# ax3.plot(t_plot, T_pc_out, linewidth = lw, label='PC HTF outlet')
# ax3.set_xlabel('Time (days)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax3.set_ylabel('Temperature (C)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax3.legend(loc=loc,fontsize=fsl)

# operating modes time history
# fig = plt.figure()
# ax = fig.gca()
# ax2.plot(t_plot, op_mode_1)
# for op in n_modes:
#     inds = op_mode_1 == op
#     ax2.plot(t_plot[inds], op_mode_1[inds], '.', label=op_modes_list[int(op)])
# ax2.set_ylabel('Operating Mode', labelpad=lp, fontsize=fs, fontweight='bold')
# ax2.legend(loc=loc,fontsize=fsl)


# fig = plt.figure(figsize=[10,4])
# ax1 = fig.gca()
# ax2 = ax1.twinx()
# # ax1 = fig.add_subplot(311)
# # ax2 = fig.add_subplot(312)
# # ax3 = fig.add_subplot(313)
# ind = 48
# rng = slice(ind*0,ind*1,1)
# t_plot = t_plot[rng]
# e_ch_tes = e_ch_tes[rng]
# p_cycle = p_cycle[rng]
# q_dot_rec_in = q_dot_rec_in[rng]
# gen = gen[rng]
# price = nuctes.df_array[rng]*200
# dt = np.diff(t_plot)[0]

# # Energy plot
# ax1.bar(t_plot,price,dt,alpha=0.5,label="Price Multiplier")
# ax2.plot(t_plot, e_ch_tes, linewidth = lw, color='C3', label='Salt Charge Level (MWh)')
# ax1.plot(t_plot, p_cycle, linewidth = lw, label='P_cycle (Electric)')
# ax1.plot(t_plot, q_dot_rec_in, linewidth = lw, label='Q_dot to Salt (Thermal)')
# ax1.plot(t_plot, gen, linewidth = lw, label='Power generated')
# ax1.set_ylabel('Power (MW)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax2.set_ylabel('Energy (MWh)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax1.legend(loc='best',fontsize=fsl)
# ax2.legend(loc=loc_cr,fontsize=fsl)
# ax1.set_xlabel('Time (days)', labelpad=lp, fontsize=fs, fontweight='bold')

# plt.tight_layout()

