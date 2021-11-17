#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:16:08 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
import os, pint
import matplotlib.pyplot as plt
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()


# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# No Dispatch Targets
# =============================================================================

# initialize the PySSC Wrapper
pw = PySSCWrapper(json_name='model1',is_debug=True)

# update SSC dictionary parameters
pw.sscdict['tshours']   = 2
pw.sscdict['P_ref']     = 450

# run SSC through PySSC
pw.run_sim(run_dispatch_targets=False)

# =============================================================================
# Extract Data
# =============================================================================
t_plot        = pw.get_array('time_hr')
last_ind = int( t_plot[t_plot>0][-1] )
out_slice = slice(0,last_ind,1)

# get output arrays from the SSC run
t_plot = t_plot[out_slice] / 24
p_cycle       = pw.get_array('P_cycle')[out_slice]
q_cycle       = pw.get_array('q_pb')[out_slice]
gen           = pw.get_array('gen')[out_slice] / 1e3
p_cool        = pw.get_array('P_cooling_tower_tot')[out_slice]
q_dot_nuc_in  = pw.get_array('q_dot_nuc_inc')[out_slice]
Q_nuc_thermal = pw.get_array('Q_nuc_thermal')[out_slice]
m_dot_nuc     = pw.get_array('m_dot_nuc')[out_slice]
T_pc_in       = pw.get_array('T_pc_in')[out_slice]
T_pc_out      = pw.get_array('T_pc_out')[out_slice]
e_ch_tes      = pw.get_array('e_ch_tes')[out_slice]
defocus       = pw.get_array('defocus')[out_slice]

op_mode_1     = pw.get_array('op_mode_1')[out_slice]

q_ch_tes  = pw.get_array('q_ch_tes')[out_slice]
q_dc_tes  = pw.get_array('q_dc_tes')[out_slice]

T_nuc_in       = pw.get_array('T_nuc_in')[out_slice]
T_nuc_out      = pw.get_array('T_nuc_out')[out_slice]

T_tes_cold_in   = pw.get_array('T_tes_cold_in')[out_slice]
T_tes_cold_out  = pw.get_array('T_tes_cold')[out_slice]
T_tes_hot       = pw.get_array('T_tes_hot')[out_slice]

q_dot_pc_target  = pw.get_array('q_dot_pc_target_on')[out_slice]
q_dot_pc_min     = pw.get_array('q_dot_pc_min')[out_slice]

qdbal       = pw.get_array('q_balance')[out_slice] * 100
mdbal       = pw.get_array('m_dot_balance')[out_slice] * 100


# =============================================================================
# Plots 
# =============================================================================

fig = plt.figure(figsize=(10,5))
ax = fig.gca()
ax.plot(t_plot, q_dot_nuc_in, '.', linewidth=2, label='qdot LFR')
ax.plot(t_plot, q_ch_tes, '.', linewidth=2, label='qdot CH')
ax.plot(t_plot, q_dc_tes, '.', linewidth=2, label='qdot DC')
ax.plot(t_plot, p_cycle, '.', linewidth=2, label='P_cycle_out')
ax.plot(t_plot, q_cycle, '.', linewidth=2, label='Q_cycle_in')
ax.plot(t_plot, q_dot_pc_target, '--', linewidth=2, label='PC thermal rating')
ax.plot(t_plot, q_dot_pc_min, '--', linewidth=2, label='PC min')
# if 
ax.set_xlabel('Time (d)', fontweight='bold')
ax.set_ylabel('Qdot (MW)', fontweight='bold')
ax.legend(loc='upper right')

# =============================================================================
# Defocus

fig = plt.figure()
ax = fig.gca()
ax.plot(t_plot, defocus, '.', linewidth=2)
ax.set_ylim([-.1,1.1])
ax.set_xlabel('Time (d)', fontweight='bold')
ax.set_ylabel('Defocus Fraction', fontweight='bold')

# =============================================================================
# Error


fig = plt.figure()
ax = fig.gca()
ax.plot(t_plot, qdbal, '.', linewidth=2, label='qdot error')
ax.plot(t_plot, mdbal, '.', linewidth=2, label='mdot error')
ax.set_xlabel('Time (d)', fontweight='bold')
ax.set_ylabel('Error Percent (%)', fontweight='bold')
ax.legend(loc='upper right')

# =============================================================================
# OP modes

# plotting OP modes
op_mode_result, modes_order = np.unique(op_mode_1, return_index=True) # mode orders and re-ordering
op_mode_result = op_mode_result[np.argsort(modes_order)]  

operating_modes = [
    'ITER_START',
    'CR_OFF__PC_OFF__TES_OFF',
    'CR_SU__PC_OFF__TES_OFF',
    'NUC_ON__PC_SU__TES_OFF',
    'NUC_ON__PC_SB__TES_OFF',
    'NUC_ON__PC_RM_HI__TES_OFF',
    'NUC_ON__PC_RM_LO__TES_OFF',
    'CR_DF__PC_MAX__TES_OFF',
    'CR_OFF__PC_SU__TES_DC',
    'NUC_ON__PC_OFF__TES_CH',
    'NUC_ON__PC_TARGET__TES_CH',
    'NUC_ON__PC_TARGET__TES_DC',
    'NUC_ON__PC_RM_LO__TES_EMPTY',
    'CR_DF__PC_OFF__TES_FULL',
    'CR_OFF__PC_SB__TES_DC',
    'CR_OFF__PC_MIN__TES_EMPTY',
    'CR_OFF__PC_RM_LO__TES_EMPTY',
    'NUC_ON__PC_SB__TES_CH',
    'CR_SU__PC_MIN__TES_EMPTY',
    'CR_SU__PC_SB__TES_DC',
    'NUC_ON__PC_SB__TES_DC',
    'CR_OFF__PC_TARGET__TES_DC',
    'CR_SU__PC_TARGET__TES_DC',
    'NUC_ON__PC_RM_HI__TES_FULL',
    'NUC_ON__PC_MIN__TES_EMPTY',
    'CR_SU__PC_RM_LO__TES_EMPTY',
    'CR_DF__PC_MAX__TES_FULL',
    'NUC_ON__PC_SB__TES_FULL',
    'CR_SU__PC_SU__TES_DC',
    'NUC_ON__PC_SU__TES_CH',
    'CR_DF__PC_SU__TES_FULL',
    'CR_DF__PC_SU__TES_OFF',
    'CR_TO_COLD__PC_TARGET__TES_DC',
    'CR_TO_COLD__PC_RM_LO__TES_EMPTY',
    'CR_TO_COLD__PC_SB__TES_DC',
    'CR_TO_COLD__PC_MIN__TES_EMPTY',
    'CR_TO_COLD__PC_OFF__TES_OFF',
    'CR_TO_COLD__PC_SU__TES_DC',
    'ITER_END']

fig = plt.figure()
ax = fig.gca()

for op in op_mode_result:
    # getting unique operating modes
    inds = (op_mode_1 == op)
    # individual index getting plotted with unique color and label
    if np.sum(inds):
        ax.plot(t_plot[inds], op_mode_1[inds], 'o', label=operating_modes[int(op)])
ax.legend()