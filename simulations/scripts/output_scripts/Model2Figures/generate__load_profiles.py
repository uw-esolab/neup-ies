#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 11:50:38 2022

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


for case in range(2,6):

    if case%2 == 0:
        # modifying inputs
        json = "model2_Hamilton_560_tariffx1_mod"
    else:
        json = "model2_CAISO_Hamilton_mod"
    dispatch = True
    run_loop = True
    sscH    = 24   # (hr)
    pyoH    = 48   # (hr)
    
    if case<=1:
        q_dot_nuclear_des=1900
    elif case<=3:
        q_dot_nuclear_des=950
    else:
        q_dot_nuclear_des=250
    
    if case==0:
        Pref    = 1100 # (MW)
        tshours = 1    # (hr)
    elif case==1:
        Pref = 1300
        tshours = 5
    elif case==2:
        Pref = 650
        tshours = 1
    elif case==3:
        Pref = 750
        tshours = 6
    elif case==4:
        Pref = 300
        tshours = 4
    else:
        Pref = 375
        tshours = 6
        
    
    
    # ========================
    
    # defining directories
    nuctes = DualPlantTES.DualPlantTES(json_name=json, is_dispatch=dispatch, log_dispatch_targets=False) 
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
    nuctes.SSC_dict['q_dot_nuclear_des']=q_dot_nuclear_des
    nuctes.dispatch_wrap.set_design()
    
    print(json )
    print("Pref    = {0}".format(nuctes.SSC_dict['P_ref'] ) )
    print("tshours = {0}".format(nuctes.SSC_dict['tshours'] ) )
    # ========================
    
    nuctes.run_sim( run_loop=run_loop, export=False, filename=output_file )
    nt = nuctes.Plant
    so = nuctes.SO
    
    print('Made it past execute.')
    
    
    # =============================================================================
    #     Display Results
    # =============================================================================
    mod_out = nt.PySAM_Outputs
    
    # storage dictionary
    Storage = {}
    Storage['time']  = np.asarray(mod_out.time_hr) * u.hr
    Storage['p_cycle']  = np.asarray(mod_out.P_cycle) * u.MW
    Storage['gen']      = (np.asarray(mod_out.gen) * u.kW).to('MW')
    Storage['q_nuc_thermal']  = np.asarray(mod_out.Q_nuc_thermal) * u.MW
    Storage['q_pb']         = np.asarray(mod_out.q_pb) * u.MW
    Storage['q_dot_pc_su']  = np.asarray(mod_out.q_dot_pc_startup) * u.MW
    Storage['price']    = np.asarray(nt.TimeOfDeliveryFactors.dispatch_factors_ts)
    Storage['e_ch_tes']  = np.asarray(mod_out.e_ch_tes) * u.MWh
    
    # locating output directory
    output_dir = os.path.join( FileMethods.output_dir, "" )
    filename = 'loadProfile_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_{4}__PC_{5}.nuctes'.format(
                    json, dispatch, sscH, pyoH, tshours, Pref )
    NTPath = os.path.join(output_dir, filename)
    
    # pickling
    with open(NTPath, 'wb') as f:
        pickle.dump(Storage, f)


