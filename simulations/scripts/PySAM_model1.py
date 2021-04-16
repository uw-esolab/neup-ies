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
nuctes.run_sim( run_loop=False )
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

p_cycle       = np.asarray( nt.Outputs.P_cycle )
gen           = np.asarray( nt.Outputs.gen ) / 1e3
p_cool        = np.asarray( nt.Outputs.P_cooling_tower_tot )
q_dot_rec_in  = np.asarray( nt.Outputs.q_dot_rec_inc )
m_dot         = np.asarray( nt.Outputs.m_dot_rec )
T_pc_in       = np.asarray( nt.Outputs.T_pc_in )
T_pc_out      = np.asarray( nt.Outputs.T_pc_out )
e_ch_tes      = np.asarray( nt.Outputs.e_ch_tes )
t_plot        = np.asarray( nt.Outputs.time_hr ) / 24
op_mode_1     = np.asarray( nt.Outputs.op_mode_1 )

#operating modes, copied from ssc/tcs/csp_solver.cpp
n_modes, modes_order = np.unique(op_mode_1,return_index=True)
n_modes = n_modes[np.argsort(modes_order)] # re-order modes by first appearance of each
op_modes_list = [
    "ENTRY_MODE",
    "CR_OFF__PC_OFF__TES_OFF,"
    "CR_SU__PC_OFF__TES_OFF",
    "CR_ON__PC_SU__TES_OFF",
    "CR_ON__PC_SB__TES_OFF",        
    "CR_ON__PC_RM_HI__TES_OFF",
    "CR_ON__PC_RM_LO__TES_OFF",        
    "CR_DF__PC_MAX__TES_OFF",
    "CR_OFF__PC_SU__TES_DC",
    "CR_ON__PC_OFF__TES_CH",
    "SKIP_10",
    "CR_ON__PC_TARGET__TES_CH",
    "CR_ON__PC_TARGET__TES_DC",
    "CR_ON__PC_RM_LO__TES_EMPTY",
    "CR_DF__PC_OFF__TES_FULL",       
    "CR_OFF__PC_SB__TES_DC",
    "CR_OFF__PC_MIN__TES_EMPTY",
    "CR_OFF__PC_RM_LO__TES_EMPTY",
    "CR_ON__PC_SB__TES_CH",
    "CR_SU__PC_MIN__TES_EMPTY",
    "SKIP_20",
    "CR_SU__PC_SB__TES_DC",
    "CR_ON__PC_SB__TES_DC",
    "CR_OFF__PC_TARGET__TES_DC",
    "CR_SU__PC_TARGET__TES_DC",
    "CR_ON__PC_RM_HI__TES_FULL",
    "CR_ON__PC_MIN__TES_EMPTY",
    "CR_SU__PC_RM_LO__TES_EMPTY",
    "CR_DF__PC_MAX__TES_FULL",
    "CR_ON__PC_SB__TES_FULL",
    "SKIP_30",
    "CR_SU__PC_SU__TES_DC",
    "CR_ON__PC_SU__TES_CH",
    "CR_DF__PC_SU__TES_FULL",
    "CR_DF__PC_SU__TES_OFF",
    "CR_TO_COLD__PC_TARGET__TES_DC",
    "CR_TO_COLD__PC_RM_LO__TES_EMPTY",
    "CR_TO_COLD__PC_SB__TES_DC",
    "CR_TO_COLD__PC_MIN__TES_EMPTY",
    "CR_TO_COLD__PC_OFF__TES_OFF",
    "SKIP_40",
    "CR_TO_COLD__PC_SU__TES_DC" ]

lp = 16 #labelpad
fs = 12 #fontsize
lw = 2  #linewidth
fsl = 'x-small'      #fontsize legend
loc = 'upper right'  #location of legend

fig = plt.figure(figsize=[10,8])
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

# Energy plot
ax1.plot(t_plot, e_ch_tes, linewidth = lw, label='Salt Charge Level (MWht)')
ax1.plot(t_plot, p_cycle, linewidth = lw, label='P_cycle (Electric)')
ax1.plot(t_plot, q_dot_rec_in, linewidth = lw, label='Q_dot to Salt (Thermal)')
ax1.plot(t_plot, gen, linewidth = lw, label='Power generated')
ax1.set_ylabel('Power (MW)', labelpad=lp, fontsize=fs, fontweight='bold')
ax1.legend(loc=loc,fontsize=fsl)

# Mass flow rate plot
# ax2.plot(t_plot, m_dot, linewidth = lw, label='m_dot_water_pc')
# ax2.set_ylabel('Mass flow (kg/s)', labelpad=lp, fontsize=fs, fontweight='bold')
# ax2.legend(loc=loc,fontsize=fsl)

# Temperature plot
ax3.plot(t_plot, T_pc_in, linewidth = lw, label='PC HTF inlet')
ax3.plot(t_plot, T_pc_out, linewidth = lw, label='PC HTF outlet')
ax3.set_xlabel('Time (days)', labelpad=lp, fontsize=fs, fontweight='bold')
ax3.set_ylabel('Temperature (C)', labelpad=lp, fontsize=fs, fontweight='bold')
ax3.legend(loc=loc,fontsize=fsl)

# operating modes time history
# fig = plt.figure()
# ax = fig.gca()
ax2.plot(t_plot, op_mode_1)
for op in n_modes:
    inds = op_mode_1 == op
    ax2.plot(t_plot[inds], op_mode_1[inds], '.', label=op_modes_list[int(op)])
ax2.set_ylabel('Operating Mode', labelpad=lp, fontsize=fs, fontweight='bold')
ax2.legend(loc=loc,fontsize=fsl)