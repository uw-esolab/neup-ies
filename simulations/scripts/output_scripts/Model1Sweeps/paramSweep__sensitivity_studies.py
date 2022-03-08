#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 15:23:51 2022

@author: gabrielsoto
"""


import modules.NuclearTES as NuclearTES
from util.PostProcessing import OutputExtraction
from util.FileMethods import FileMethods
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import numpy as np
import hashlib
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# Static  - User Defined Tables + Tariff Price Schedule
# Varying - PC size, TES size
# =============================================================================

# =============================================================================
# Parameters for Sweep
# =============================================================================
sscH = 24  # 12 # 24
pyoH = 48  # 24 # 48
json = "model1_Hamilton_560_tariffx1_5"  # model1_CAISO_Hamilton # model1 # model1_Hamilton_560_tariffx2 
sensitivityStudy = "TESCost" # HotTemp # TESCost
dispatch = True # True # False
run_loop = True

# costs and other params
turb_unit_cost = 225 #Cory email 8th Nov 2021
nuc_spec_cost  = 4150
tes_spec_cost  = 28.4
fin_yrs  = 4.0
fin_rate = 0.07

tes_costs = np.array([   20,  30,  40,  50,  60])
hot_temp  = np.array([  550, 555, 560, 565, 570])

if json is "model1_CAISO_Hamilton":
    tshours   = np.array([    2,   4,   6 ])
    p_cycle   = np.array([ 1100, 950, 800 ]) 
    
elif json is "model1_Hamilton_560_tariffx1":
    tshours   = np.array([    0,   1,   2 ])
    p_cycle   = np.array([  550, 500, 450 ]) 
    
elif json is "model1_Hamilton_560_tariffx1_5" or "model1_Hamilton_560_tariffx2":
    tshours   = np.array([    4,   6,   8 ])
    p_cycle   = np.array([  950, 800, 650 ]) 

if sensitivityStudy is "TESCost":
    iterator0 = tes_costs
elif sensitivityStudy is "HotTemp":
    iterator0 = hot_temp


# =============================================================================
# Initialize empty arrays
# =============================================================================

# formally defining the iterator arrays
iterator1 = tshours  
iterator2 = p_cycle

# initializing an empty array
empty   = np.zeros([  len(iterator0) , len(iterator1) , len(iterator2) ])
empty2D = np.zeros([  len(iterator0) , len(iterator1) , len(iterator2), 2 ])

# initializing output arrays
defocus   = empty2D.copy()
gen       = empty2D.copy()
qdot_nuc  = empty2D.copy()
revenue   = empty.copy()
ppa       = empty.copy()
annual_e  = empty.copy()
cap_fac   = empty.copy()
nuc_cost  = empty.copy()
fin_cost  = empty.copy()
TES_CH    = empty.copy()
TES_DC    = empty.copy()
iter_log  = empty.copy()
fail_log  = empty.copy() # 0 = No failure, 1 = SSC caused failure, 
                         # 2 = Pyomo caused failure, 3 = Pyomo crash occurred
                         # 4 = SSC low mass flow error
pyomo_bad_log     = empty.copy()
pyomo_bad_idx_log = empty.copy()
pysamdict = {}
# =============================================================================
# Operating Modes analysis
# =============================================================================

from util.PostProcessing import Plots
tmp_modes = lambda: None
Plots.set_operating_modes_list(tmp_modes)

op_modes_list = tmp_modes.operating_modes


# =============================================================================
# Double Loop
# =============================================================================

# starting the time counter
tic = time.perf_counter()

# TRIPLE LOOP
for s, sp in enumerate(iterator0): #over sensitivity parameter
    for i, th in enumerate(iterator1): #over tshours
        for j, pc in enumerate(iterator2): #over P_cycle rating
            
            # =========================================
            # set up
            # =========================================
            
            # getting indeces
            idx = (s,i,j)
            mean_idx = (s,i,j,0)
            stdv_idx = (s,i,j,1)
            
            # =========================================
            # updating pysam class
            # =========================================
            
            # defining directories
            nuctes = NuclearTES.NuclearTES( json_name=json, is_dispatch=dispatch)
    
            # horizons
            nuctes.PySAM_dict['ssc_horizon']   = sscH
            nuctes.ssc_horizon   = sscH * nuctes.u.hr
            nuctes.PySAM_dict['pyomo_horizon'] = pyoH
            nuctes.pyomo_horizon = pyoH * nuctes.u.hr
            
            # saving/updating PYSAM dict to nuctes
            nuctes.dispatch_wrap = nuctes.create_dispatch_wrapper( nuctes.PySAM_dict )
    
            # update SSC dictionary parameters
            nuctes.SSC_dict['P_ref']   = float(pc)
            nuctes.SSC_dict['tshours'] = float(th)
            
            # update turbine costs
            turb_ref = nuctes.SSC_dict['q_dot_nuclear_des'] * nuctes.SSC_dict['design_eff']
            turb_premium = turb_unit_cost * (nuctes.SSC_dict['P_ref'] - turb_ref)  
            nuctes.SSC_dict["nuclear_spec_cost"] = (nuc_spec_cost*turb_ref + turb_premium) / nuctes.SSC_dict['P_ref']
            
            # update TES costs
            if sensitivityStudy is "TESCost":
                nuctes.SSC_dict["tes_spec_cost"] = float(sp)
            else:
                nuctes.SSC_dict["tes_spec_cost"] = tes_spec_cost
            
            if sensitivityStudy is "HotTemp":
                nuctes.SSC_dict["T_htf_hot_des"]   = float(sp)
                nuctes.SSC_dict["T_tank_hot_init"] = float(sp)
            
            
            # update construction financial costs
            base_financing = nuctes.SSC_dict["construction_financing_cost"]
            extra_financing = fin_yrs * fin_rate * 1000 * (turb_premium + \
                                                    nuctes.SSC_dict['P_ref'] * nuctes.SSC_dict['tshours'] * nuctes.SSC_dict["tes_spec_cost"])
            nuctes.SSC_dict["construction_financing_cost"]=base_financing + extra_financing
                        
            # reset design point in Dispatch parameter class
            nuctes.dispatch_wrap.set_design()
            
            # =========================================
            # actually run this simulation
            # =========================================
            
            # print current position in loop
            print(json)
            print("TEScost = {0}".format(nuctes.SSC_dict['tes_spec_cost'] ) )
            print("HotTemp = {0}".format(nuctes.SSC_dict['T_htf_hot_des'] ) )
            print("Pref    = {0}".format(nuctes.SSC_dict['P_ref'] ) )
            print("tshours = {0}".format(nuctes.SSC_dict['tshours'] ) )
            
            # =========================================
            # try running simulation, log success
            # =========================================
            
            success = True
            try:
                # run SSC through PySSC
                nuctes.run_sim( run_loop=run_loop )
            except Exception as err:
                print(" Run failed.")
                success = False
    
            # =========================================
            # log success of simulation
            # =========================================
            
            # extract Plant object
            nt = nuctes.Plant
            
            # log final time index of simulation
            iter_log[idx] = nuctes.slice_ssc_currentH.start if nuctes.slice_ssc_currentH.stop != 8760 \
                                    else nuctes.slice_ssc_currentH.stop
            
            # log success
            success = False if nuctes.slice_ssc_currentH.stop != 8760 else success
            
            # if we ran dispatch optimization, check for failures or crashes 
            if dispatch:
                
                # log all of the Pyomo failures (the dictionary stores successes, we negate it with the 'not' command)
                pyomo_fail = np.array([not nuctes.disp_success[k] for k in nuctes.disp_success.keys()])
                
                # number of failures due to Pyomo
                pyomo_fail_sum = sum(pyomo_fail)
                
                # log whether Pyomo failed at all throughout simulation
                pyomo_bad_log[idx] = bool(pyomo_fail_sum > 0)
                
                # log the time index where Pyomo crashed if it did indeed crash
                if pyomo_fail_sum > 0:
                    pyomo_bad_idx_log[idx] = np.where(pyomo_fail==True)[0][0]
    
            # =========================================
            # Simulation was successful! Log results
            # =========================================
            if success:
                so = nuctes.SO
                outputs = OutputExtraction(nuctes)
                
                # analyzing OP modes
                op_mode_array       = outputs.op_mode_1
                op_mode_str_profile = [op_modes_list[ int(x) ] for x in op_mode_array]
                
                TES_charging    = np.array(["TES_CH" in x for x in op_mode_str_profile])
                TES_discharging = np.array(["TES_DC" in x for x in op_mode_str_profile])
                
                TES_CH[idx] = int(TES_charging.sum() > 0)
                TES_DC[idx] = int(TES_discharging.sum() > 0)
                
                # log output means
                gen[mean_idx]       = np.mean( outputs.gen.m )
                qdot_nuc[mean_idx]  = np.mean( outputs.q_dot_nuc_in.m )
                defocus[mean_idx]   = np.mean( outputs.defocus )
                
                # log output std dev
                gen[stdv_idx]       = np.std( outputs.gen.m )
                qdot_nuc[stdv_idx]  = np.std( outputs.q_dot_nuc_in.m )
                defocus[stdv_idx]   = np.std( outputs.defocus )
                
                # log economic outputs
                revenue[idx]   = np.sum( outputs.gen.to('kW').m * outputs.price )
                ppa[idx]       = so.Outputs.ppa
                annual_e[idx]  = (nt.Outputs.annual_energy*u.kWh).to('TWh').m # TWh
                cap_fac[idx]   = nuctes.capacity_factor.m * 100
                nuc_cost[idx]  = nt.SystemCosts.nuclear_spec_cost
                fin_cost[idx]  = so.FinancialParameters.construction_financing_cost
                
                print("ppa  = {0} cents/kWh \n\n".format( ppa[idx] ) )
                print("cost = {0} cents/kWh \n\n".format( fin_cost[idx] ) )
                
                del outputs, op_mode_array, op_mode_str_profile, so
            
            # =========================================
            # Simulation failed! Log results
            # =========================================
            else:
                
                # check if last crash was due to Pyomo
                if dispatch:
                    # checking the last success log at end of simulation
                    count = int(nuctes.disp_count - 1)
                    disp_success = nuctes.disp_success[count]
                    
                    # Simulation ended because of SSC crash   == 1 
                    # Simulation ended because of Pyomo crash == 2 
                    fail_log[idx] = 1 if disp_success else 2 
                
                if hasattr(nuctes,'err_message'):
                    if 'Mass flow rate too low' in nuctes.err_message:
                        fail_log[idx] = 4
                    elif 'stuck' in nuctes.err_message:
                        fail_log[idx] = 5
                    
                
                dummy_val = -1
                # log output means
                gen[mean_idx]       = dummy_val
                qdot_nuc[mean_idx]  = dummy_val
                defocus[mean_idx]   = dummy_val
                
                # log output std dev
                gen[stdv_idx]       = dummy_val
                qdot_nuc[stdv_idx]  = dummy_val
                defocus[stdv_idx]   = dummy_val
                
                # log economic results
                revenue[idx]   = dummy_val
                ppa[idx]       = dummy_val
                annual_e[idx]  = dummy_val
                cap_fac[idx]   = dummy_val
                nuc_cost[idx]  = dummy_val
                fin_cost[idx]  = dummy_val
                
                TES_CH[idx] = dummy_val
                TES_DC[idx] = dummy_val 
                    
                # reset the Plant and Grid, prevents memory leak
            del nt
            
            pysamdict = nuctes.PySAM_dict
            
            del nuctes
                
# end time counter
toc = time.perf_counter()        

# print out time elapsed and store it (current record is 5.6 hours)
time_elapsed = ((toc-tic)*u.s).to('hr')
print('Made it past the triple loop. Time %.2f hrs' % time_elapsed.m)

# =============================================================================
# Storing Data
# =============================================================================

# storage dictionary
Storage = {}
Storage['time_elapsed']  = time_elapsed
Storage['tshours']    = tshours
Storage['p_cycle']    = p_cycle
Storage['sensitivityStudy'] = sensitivityStudy
Storage['iterator0']  = iterator0
Storage['iterator1']  = iterator1
Storage['iterator2']  = iterator2
Storage['sscH']       = sscH
Storage['pyoH']       = pyoH
Storage['op_modes_list'] = op_modes_list
Storage['json']       = json
Storage['dispatch']   = dispatch
Storage['run_loop']   = run_loop
Storage['defocus']    = defocus
Storage['gen']        = gen
Storage['qdot_nuc']   = qdot_nuc
Storage['TES_CH']     = TES_CH
Storage['TES_DC']     = TES_DC
Storage['iter_log']   = iter_log
Storage['fail_log']   = fail_log
Storage['pyomo_bad_log']      = pyomo_bad_log
Storage['pyomo_bad_idx_log']  = pyomo_bad_idx_log
Storage['revenue']  = revenue
Storage['annual_e'] = annual_e
Storage['ppa']      = ppa
Storage['cap_fac']  = cap_fac
Storage['turb_unit_cost']    = turb_unit_cost
Storage['tes_spec_cost']     = tes_spec_cost
Storage['fin_yrs']     = fin_yrs
Storage['fin_rate']    = fin_rate
Storage['nuclear_spec_cost']            = nuc_cost
Storage['construction_financing_cost'] = fin_cost


coeff_list = ['ec_p',     # proft term
              'ec_pcsu',  # PC cold SU term
              'ec_pchsu', # PC hot SU term 
              'ec_pcsd',  # PC SD term
              'ec_pcr',   # PC ramping term
              'ec_pcer',  # PC excess ramping term
              'ec_nsu',   # LFR cold SU term
              'ec_nhsu',  # LFR hot SU term
              'ec_nsd',   # LFR SD term
              'ec_pcg',   # PC power generating term
              'ec_pcsb',  # PC SB term
              'ec_nt'     # LFR thermal power generating term
              ]
        
# setting objective function coefficients
# extr_str = ''
# for coeff in coeff_list:
#     Storage[coeff] = pysamdict[coeff] if coeff in pysamdict.keys() else 1.0
#     extr_str += "_{0}{1}".format( coeff.split('_')[1], str(Storage[coeff]))

# locating output directory
output_dir = FileMethods.output_dir
filename = 'paramSweep_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_[{4},{5}]__PC_[{6},{7}]__sens_{8}.nuctes'.format(
                json, dispatch, sscH, pyoH, tshours.min(), tshours.max(), p_cycle.min(), p_cycle.max(), sensitivityStudy )
NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'wb') as f:
    pickle.dump(Storage, f)