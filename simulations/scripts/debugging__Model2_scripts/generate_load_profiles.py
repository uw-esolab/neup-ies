#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:29:24 2022

@author: gabrielsoto
"""


import os,sys
sys.path.append('..')
import modules.DualPlantTES as DualPlantTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, time, copy, pickle
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
json = "model2_Hamilton_560_TwinPeaks_x1"   # model2_CAISO_Hamilton # model2_Hamilton_560_tariffx1 # model2_Hamilton_560_TwinPeaks_x1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 800 # (MW)
tshours = 8    # (hr)
qdotrec = 500  # (MW)

# ========================

# defining NE2 module
dptes = DualPlantTES.DualPlantTES(json_name=json, is_dispatch=dispatch) 

# horizons
dptes.PySAM_dict['ssc_horizon']   = sscH
dptes.ssc_horizon   = sscH * dptes.u.hr
dptes.PySAM_dict['pyomo_horizon'] = pyoH
dptes.pyomo_horizon = pyoH * dptes.u.hr

# saving/updating PYSAM dict to nuctes
dptes.dispatch_wrap = dptes.create_dispatch_wrapper( dptes.PySAM_dict )

# update SSC dictionary parameters
dptes.SSC_dict['P_ref'] = Pref
dptes.SSC_dict['tshours'] = tshours
dptes.SSC_dict['q_dot_rec_des'] = qdotrec
# dptes.SSC_dict['q_sby_frac'] = 1e-5
# dptes.SSC_dict['cycle_cutoff_frac'] = 1e-5
# dptes.SSC_dict['rec_qf_delay'] = 1e-5
# dptes.SSC_dict['f_rec_min'] = 1e-5

dptes.dispatch_wrap.set_design()

print("Pref    = {0}".format(dptes.SSC_dict['P_ref'] ) )
print("tshours = {0}".format(dptes.SSC_dict['tshours'] ) )
print("QdotRec = {0}".format(dptes.SSC_dict['q_dot_rec_des'] ) )
# ========================

dptes.run_sim( run_loop=run_loop, export=False )
dt = dptes.Plant
so = dptes.SO

print('Made it past execute.')


# =============================================================================
#     Display Results
# =============================================================================
mod_out = dt.PySAM_Outputs

# storage dictionary
Storage = {}
Storage['time']  = np.asarray(mod_out.time_hr) * u.hr
Storage['p_cycle']  = np.asarray(mod_out.P_cycle) * u.MW
Storage['gen']      = (np.asarray(mod_out.gen) * u.kW).to('MW')
Storage['q_nuc_thermal']  = np.asarray(mod_out.Q_nuc_thermal) * u.MW
Storage['q_rec_thermal']  = np.asarray(mod_out.Q_thermal) * u.MW
Storage['q_pb']         = np.asarray(mod_out.q_pb) * u.MW
Storage['q_dot_pc_su']  = np.asarray(mod_out.q_dot_pc_startup) * u.MW
Storage['price']    = np.asarray(dt.TimeOfDeliveryFactors.dispatch_factors_ts)
Storage['e_ch_tes']  = np.asarray(mod_out.e_ch_tes) * u.MWh

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "model2a" )
filename = 'loadProfile_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_{4}__PC_{5}__CSP_{6}.dualtes'.format(
                json, dispatch, sscH, pyoH, tshours, Pref, qdotrec )
NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'wb') as f:
    pickle.dump(Storage, f)