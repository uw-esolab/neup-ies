#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:27:13 2022

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
import os

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# initialize the PySSC Wrapper
pw = PySSCWrapper(json_name='model2b_Hamilton_FS_tariffx1',is_debug=True)

# update SSC dictionary parameters
pw.sscdict['P_ref']         = 700
pw.sscdict['tshours']       = 8
pw.sscdict['q_dot_rec_des'] = 500

# run SSC through PySSC
pw.run_sim(run_dispatch_targets=True)

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
