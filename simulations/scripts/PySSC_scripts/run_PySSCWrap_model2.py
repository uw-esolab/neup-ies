#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 09:56:04 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
import os

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# initialize the PySSC Wrapper
pw = PySSCWrapper(json_name='model2',is_debug=True)

# run SSC through PySSC
pw.run_sim()

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
