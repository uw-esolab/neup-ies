#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 13:14:51 2022

@author: gabrielsoto
"""


import modules.DualPlantTES as DualPlantTES
from util.DualPlots import DualOutputExtraction
from util.FileMethods import FileMethods
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# set up
# =============================================================================

# possible different CSP design qdot 
csp_qdot_ref = 750 # 500 # 400

# array of successful P_cycle and tshours combinations
# p_cycle_0     = np.array([ 700, 650, 600, 550,      450, 400, 350, 300, 250, 200, 150 ]) # Pref= 500, tshours=0
# p_cycle_2     = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=2
# p_cycle_4     = np.array([ 700, 650, 600, 550, 500,                300, 250, 200, 150 ]) # Pref= [350,450], tshours=4
# p_cycle_6     = np.array([ 700, 650, 600, 550, 500,      400,      300, 250, 200, 150 ]) # Pref= {350,450}, tshours=6
# p_cycle_8     = np.array([ 700, 650, 600, 550, 500,           350, 300, 250,      150 ]) # Pref= {200,400,450}, tshours=8
# p_cycle_10    = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=10
# p_cycle_12    = np.array([ 700, 650, 600, 550, 500,           350, 300, 250,      150 ]) # Pref= {200,400,450}, tshours=12
# p_cycle_14    = np.array([ 700, 650, 600, 550, 500,      400, 350, 300, 250, 200, 150 ]) # Pref= 450, tshours=14

# setting up arrays to cycle through
# tshours    = np.array([ 0, 2, 4, 6, 8, 10, 12 ])
tshours    = np.array([ 0, 2, 4, 6, 8, 10, 12, 14 ])
# p_cycle    = np.array([ 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100 ]) 
# p_cycle    = np.array([ 500, 450, 400, 350, 300, 250, 200, 150, 100 ]) 
# p_cycle    = np.array([ 500, 450, 400, 350, 300, 250, 200 ]) 
# p_cycle    = np.array([ 850, 800, 750, 700, 650, 600, 550 ]) 

# p_cycle    = np.array([ 1100, 1050, 1000, 950, 900]) 
p_cycle    = np.array([ 850, 800, 750, 700, 650, 600]) 
# p_cycle    = np.array([ 550, 500, 450, 400, 350, 300 ]) 
# 
# p_cycle    = np.array([ 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300 ]) 

# exceptions!
# exceptions = {
#     "0":  [500],
#     "2":  [450],
#     "4":  [350, 400, 450],
#     "6":  [350, 450],
#     "8":  [200, 400, 450],
#     "10": [450],
#     "12": [200, 400, 450], 
#     "14": [450]
#     }

# exceptions = {"0":  [500, 450, 400], "2":  [400], "4":  [], "6":  [400], "8":  [], "10": [400], "12": [], "14": [400] }
exceptions = {"0":  [], "2":  [], "4":  [], "6":  [], "8":  [], "10": [], "12": [], "14": [] }

# formally defining the iterator arrays
iterator1 = tshours  
iterator2 = p_cycle

# initializing an empty array
empty   = np.zeros([  len(iterator1) , len(iterator2) ])
empty2D = np.zeros([  len(iterator1) , len(iterator2), 2 ])

# initializing output arrays
defocus   = empty2D.copy()
gen       = empty2D.copy()
qdot_nuc  = empty2D.copy()
revenue   = empty.copy()
TES_CH    = empty.copy()
TES_DC    = empty.copy()
iter_log  = empty.copy()
fail_log  = empty.copy() # 0 = No failure, 1 = SSC caused failure, 
                         # 2 = Pyomo caused failure, 3 = Pyomo crash occurred
                         # 4 = SSC low mass flow error, 5 = 
pyomo_bad_log     = empty.copy()
pyomo_bad_idx_log = empty.copy()

allOff_log  = empty.copy()

# =============================================================================
# Operating Modes analysis
# =============================================================================

from util.DualPlots import DualPlots
tmp_modes = lambda: None
DualPlots.set_operating_modes_list(tmp_modes)

op_modes_list = tmp_modes.operating_modes


# =============================================================================
# Double Loop
# =============================================================================
sscH = 24  # 12 # 24
pyoH = 48  # 24 # 48
json = "model2_Hamilton_560_TwinPeaks_x1" # model2_Hamilton_560_tariffx1 # model2_Hamilton_560_TwinPeaks_x1 # model2_CAISO_Hamilton
dispatch = True # True # False
run_loop = True


# starting the time counter
tic = time.perf_counter()

# TRIPLE LOOP
for i, th in enumerate(iterator1): #over tshours
    for j, pc in enumerate(iterator2): #over P_cycle rating
        
        # =========================================
        # set up
        # =========================================
        
        # getting indeces
        idx = (i,j)
        mean_idx = (i,j,0)
        stdv_idx = (i,j,1)

        # print current position in loop
        print(json)
        print("output tshours :   ", th)
        print("output Pref    :   ", pc)
        print("output Pref    :   ", csp_qdot_ref)
        
        # =========================================
        # updating pysam class
        # =========================================
        
        # defining directories
        dltes = DualPlantTES.DualPlantTES( json_name=json, is_dispatch=dispatch)

        # horizons
        dltes.PySAM_dict['ssc_horizon']   = sscH
        dltes.ssc_horizon   = sscH * dltes.u.hr
        dltes.PySAM_dict['pyomo_horizon'] = pyoH
        dltes.pyomo_horizon = pyoH * dltes.u.hr
        
        # saving/updating PYSAM dict to dltes
        dltes.dispatch_wrap = dltes.create_dispatch_wrapper( dltes.PySAM_dict )

        # update SSC dictionary parameters
        dltes.SSC_dict['P_ref']   = float(pc)
        dltes.SSC_dict['tshours'] = float(th)
        dltes.SSC_dict['q_dot_rec_des'] = csp_qdot_ref
        dltes.dispatch_wrap.set_design()
        
        # =========================================
        # exceptions we know don't work in SSC
        # =========================================
        
        if pc in exceptions[str(th)]:
            dummy_val = -1
            # log output means
            gen[mean_idx]       = dummy_val
            qdot_nuc[mean_idx]  = dummy_val
            defocus[mean_idx]   = dummy_val
            
            # log output std dev
            gen[stdv_idx]       = dummy_val
            qdot_nuc[stdv_idx]  = dummy_val
            defocus[stdv_idx]   = dummy_val
            
            TES_CH[idx] = dummy_val
            TES_DC[idx] = dummy_val
            
        # =========================================
        # actually run this simulation
        # =========================================

        else:
            
            # =========================================
            # try running simulation, log success
            # =========================================
            success = True
            try:
                # run SSC through PySSC
                dltes.run_sim( run_loop=run_loop )
            except Exception as err:
                print(" Run failed.")
                success = False

            # =========================================
            # log success of simulation
            # =========================================
            
            # extract Plant object
            nt = dltes.Plant
            
            # log final time index of simulation
            iter_log[idx] = dltes.slice_ssc_currentH.start if dltes.slice_ssc_currentH.stop != 8760 \
                                    else dltes.slice_ssc_currentH.stop
            
            # log success
            success = False if dltes.slice_ssc_currentH.stop != 8760 else success
            
            # if we ran dispatch optimization, check for failures or crashes 
            if dispatch:
                
                # log all of the Pyomo failures (the dictionary stores successes, we negate it with the 'not' command)
                pyomo_fail = np.array([not dltes.disp_success[k] for k in dltes.disp_success.keys()])
                
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

                outputs = DualOutputExtraction(dltes)
                
                # analyzing OP modes
                op_mode_array       = outputs.op_mode_1
                op_mode_str_profile = [op_modes_list[ int(x) ] for x in op_mode_array]
                
                # catching any OFF modes that 
                allOFF = np.array(['CR_OFF__PC_OFF__TES_OFF__NUC_ON' == x for x in op_mode_str_profile])
                allOff_log[idx] = allOFF.sum()
                
                # finding the first instance of an OFF mode that we didn't catch
                if allOFF.sum() > 0:
                    iter_log[idx] = np.where( allOFF == True )[0][0]
                
                print('All OFF: {0} \n'.format(allOff_log[idx])  )
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
                
                revenue[idx] = np.sum(outputs.gen.to('kW').m * outputs.price )
                
                del outputs, op_mode_array, op_mode_str_profile, allOFF
            
            # =========================================
            # Simulation failed! Log results
            # =========================================
            else:
                
                # check if last crash was due to Pyomo
                if dispatch:
                    # checking the last success log at end of simulation
                    count = int(dltes.disp_count - 1)
                    disp_success = dltes.disp_success[count]
                    
                    # Simulation ended because of SSC crash   == 1 
                    # Simulation ended because of Pyomo crash == 2 
                    fail_log[idx] = 1 if disp_success else 2 
                
                if hasattr(dltes,'err_message'):
                    if 'Mass flow rate too low' in dltes.err_message:
                        fail_log[idx] = 4
                    elif 'stuck' in dltes.err_message:
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
                
                TES_CH[idx] = dummy_val
                TES_DC[idx] = dummy_val 
                
            # reset the Plant and Grid, prevents memory leak
            del nt
            
        del dltes
            
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
Storage['q_dot_rec_des'] = csp_qdot_ref
Storage['exceptions'] = exceptions
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
Storage['allOff_log'] = allOff_log
Storage['pyomo_bad_log']      = pyomo_bad_log
Storage['pyomo_bad_idx_log']  = pyomo_bad_idx_log
Storage['revenue'] = revenue

# locating output directory
output_dir = FileMethods.output_dir
filename = 'failureModes_PySAM__{0}__2022_04__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_[{4},{5}]__PC_[{6},{7}]__CSP_{8}.nuctes'.format(
                json, dispatch, sscH, pyoH, tshours.min(), tshours.max(), p_cycle.min(), p_cycle.max(), csp_qdot_ref )

NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'wb') as f:
    pickle.dump(Storage, f)