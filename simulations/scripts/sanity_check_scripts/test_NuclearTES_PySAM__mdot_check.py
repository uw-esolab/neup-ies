#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 15:48:10 2021

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
Pref    = 750  # (MW)
tshours = 10    # (hr)

# ========================

# defining directories
nuctes = NuclearTES.NuclearTES(json_name=json, is_dispatch=dispatch) 
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

print("Pref    = {0}".format(nuctes.SSC_dict['P_ref'] ) )
print("tshours = {0}".format(nuctes.SSC_dict['tshours'] ) )
# ========================

nuctes.run_sim( run_loop=run_loop, export=False, filename=output_file )
nt = nuctes.Plant
so = nuctes.SO

print('Made it past execute.')


# =============================================================================
#     Plotting
# =============================================================================

from util.PostProcessing import Plots
upl = Plots(nuctes, legend_offset = True, x_shrink=0.7, fE_min=0, fE_max=1.05)


# =============================================================================
# 
# =============================================================================

from util.SSCHelperMethods import SSCHelperMethods

rec_htf = nuctes.SSC_dict['rec_htf']
T_cold  = (nuctes.T_rec_in_log * upl.u.degC).to('K')
T_hot   = (nuctes.T_rec_out_log * upl.u.degC).to('K')
T_avg   = 0.5*(T_cold+T_hot)
c_p = SSCHelperMethods.get_cp_htf( upl.u, T_avg, rec_htf )
m_dot = upl.m_dot_rec

q_calc = c_p * m_dot * (T_hot - T_cold)
q_calc = q_calc.to('MW')

fig = plt.figure()
ax = fig.gca()
ax.plot(q_calc - upl.q_dot_rec_in)

