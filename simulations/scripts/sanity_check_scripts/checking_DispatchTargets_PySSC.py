#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 15:08:17 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
import os
import matplotlib.pyplot as plt
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)


# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# No Dispatch Targets
# =============================================================================

# initialize the PySSC Wrapper
pw = PySSCWrapper(json_name='model1',is_debug=True)

# run SSC through PySSC
pw.run_sim(run_dispatch_targets=False)

# get output arrays from the SSC run
p_cycle       = pw.get_array('P_cycle')
gen           = pw.get_array('gen') / 1e3
p_cool        = pw.get_array('P_cooling_tower_tot')
q_dot_rec_in  = pw.get_array('q_dot_rec_inc')
m_dot         = pw.get_array('m_dot_rec')
T_pc_in       = pw.get_array('T_pc_in')
T_pc_out      = pw.get_array('T_pc_out')
e_ch_tes      = pw.get_array('e_ch_tes')
t_plot        = pw.get_array('time_hr') / 24
op_mode_1     = pw.get_array('op_mode_1')

# =============================================================================
# Dispatch Targets
# =============================================================================

# initialize the PySSC Wrapper
pwDT = PySSCWrapper(json_name='model1',is_debug=True)

# run SSC through PySSC
pwDT.run_sim(run_dispatch_targets=True)

# get output arrays from the SSC run
p_cycleDT     = pwDT.get_array('P_cycle')
gen           = pwDT.get_array('gen') / 1e3
p_cool        = pwDT.get_array('P_cooling_tower_tot')
q_dot_rec_in  = pwDT.get_array('q_dot_rec_inc')
m_dot         = pwDT.get_array('m_dot_rec')
T_pc_in       = pwDT.get_array('T_pc_in')
T_pc_out      = pwDT.get_array('T_pc_out')
e_ch_tes      = pwDT.get_array('e_ch_tes')
t_plot        = pwDT.get_array('time_hr') / 24
op_mode_1DT   = pwDT.get_array('op_mode_1')

fig = plt.figure()
ax = fig.gca()
ax.plot(t_plot, p_cycleDT)
ax.plot(t_plot, p_cycle)


fig = plt.figure()
ax = fig.gca()
ax.plot(t_plot, op_mode_1DT)
ax.plot(t_plot, op_mode_1)

