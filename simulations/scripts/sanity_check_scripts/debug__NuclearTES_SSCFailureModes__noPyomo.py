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
pw.sscdict['tshours']   = 0
pw.sscdict['P_ref']     = 500
# pw.sscdict['time_stop'] = (3421*u.hr).to('s').m
        
# run SSC through PySSC
pw.run_sim(run_dispatch_targets=False)

# =============================================================================
# Extract Data
# =============================================================================

# get output arrays from the SSC run
p_cycle       = pw.get_array('P_cycle')
gen           = pw.get_array('gen') / 1e3
p_cool        = pw.get_array('P_cooling_tower_tot')
q_dot_nuc_in  = pw.get_array('q_dot_nuc_inc')
m_dot_nuc     = pw.get_array('m_dot_nuc')
T_pc_in       = pw.get_array('T_pc_in')
T_pc_out      = pw.get_array('T_pc_out')
e_ch_tes      = pw.get_array('e_ch_tes')
defocus       = pw.get_array('defocus')
t_plot        = pw.get_array('time_hr') / 24
op_mode_1     = pw.get_array('op_mode_1')

q_ch_tes  = pw.get_array('q_ch_tes')
q_dc_tes  = pw.get_array('q_dc_tes')

T_nuc_in       = pw.get_array('T_nuc_in')
T_nuc_out      = pw.get_array('T_nuc_out')

T_tes_cold_in   = pw.get_array('T_tes_cold_in')
T_tes_cold_out  = pw.get_array('T_tes_cold')
T_tes_hot       = pw.get_array('T_tes_hot')


# =============================================================================
# Plots 
# =============================================================================

fig = plt.figure()
ax = fig.gca()
ax.plot(t_plot, q_dot_nuc_in, '.', linewidth=2, label='qdot LFR')
ax.plot(t_plot, p_cycle, '.', linewidth=2, label='P_cycle')
ax.plot(t_plot, gen, '.', linewidth=2, label='gen')
ax.set_xlabel('Time (d)', fontweight='bold')
ax.set_ylabel('Qdot (MW)', fontweight='bold')
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
    'CR_ON__PC_SU__TES_OFF',
    'CR_ON__PC_SB__TES_OFF',
    'CR_ON__PC_RM_HI__TES_OFF',
    'CR_ON__PC_RM_LO__TES_OFF',
    'CR_DF__PC_MAX__TES_OFF',
    'CR_OFF__PC_SU__TES_DC',
    'CR_ON__PC_OFF__TES_CH',
    'CR_ON__PC_TARGET__TES_CH',
    'CR_ON__PC_TARGET__TES_DC',
    'CR_ON__PC_RM_LO__TES_EMPTY',
    'CR_DF__PC_OFF__TES_FULL',
    'CR_OFF__PC_SB__TES_DC',
    'CR_OFF__PC_MIN__TES_EMPTY',
    'CR_OFF__PC_RM_LO__TES_EMPTY',
    'CR_ON__PC_SB__TES_CH',
    'CR_SU__PC_MIN__TES_EMPTY',
    'CR_SU__PC_SB__TES_DC',
    'CR_ON__PC_SB__TES_DC',
    'CR_OFF__PC_TARGET__TES_DC',
    'CR_SU__PC_TARGET__TES_DC',
    'CR_ON__PC_RM_HI__TES_FULL',
    'CR_ON__PC_MIN__TES_EMPTY',
    'CR_SU__PC_RM_LO__TES_EMPTY',
    'CR_DF__PC_MAX__TES_FULL',
    'CR_ON__PC_SB__TES_FULL',
    'CR_SU__PC_SU__TES_DC',
    'CR_ON__PC_SU__TES_CH',
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