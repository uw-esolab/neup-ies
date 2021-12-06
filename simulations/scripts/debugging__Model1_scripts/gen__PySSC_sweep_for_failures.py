#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 22:17:06 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
from util.PostProcessing import OutputExtraction
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
tshours    = np.array([ 0, 2, 4, 6, 8, 10, 12, 14 ])
p_cycle    = np.array([ 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100 ]) 

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
TES_CH    = empty.copy()
TES_DC    = empty.copy()
iter_log  = empty.copy()
fail_log  = empty.copy() # 0 = No failure, 1 = SSC caused failure, 2 = Pyomo caused failure, 3 = Pyomo crash occurred, but sim continued
pyomo_bad_log     = empty.copy()
pyomo_bad_idx_log = empty.copy()

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
sscH = 24  # 12 # 24
pyoH = 48  # 24 # 48
json = "model1_Hamilton_560_dfe" # model1_CAISO # model1 # model1_noMin # model1_Hamilton_560
dispatch = False # True # False
run_loop = False


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
        print("output tshours :   ", th)
        print("output Pref    :   ", pc)
        
        # =========================================
        # updating pysam class
        # =========================================
        
        # defining directories
        pw = PySSCWrapper(json_name='model1',is_debug=False)

        # horizons
        pw.sscdict['tshours'] = float(th)
        pw.sscdict['P_ref']   = float(pc)
        
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
                pw.run_sim()
            except Exception as err:
                print(" Run failed.")
                success = False

            # =========================================
            # log success of simulation
            # =========================================
            t_plot = pw.get_array('time_hr')
            
            # log final time index of simulation
            last_ind = int( t_plot[t_plot>0][-1] )
            out_slice = slice(0,last_ind,1)
            
            print( "      last index: {0}".format(last_ind) )
            
            iter_log[idx] = last_ind

            # log success
            success = False if last_ind != 8760 else success
            
            # if we ran dispatch optimization, check for failures or crashes 
            if dispatch:
                pass

            # =========================================
            # Simulation was successful! Log results
            # =========================================
            if success:

                # analyzing OP modes
                op_mode_array       = pw.get_array('op_mode_1')[out_slice]
                op_mode_str_profile = [op_modes_list[ int(x) ] for x in op_mode_array]
                
                TES_charging    = np.array(["TES_CH" in x for x in op_mode_str_profile])
                TES_discharging = np.array(["TES_DC" in x for x in op_mode_str_profile])
                
                TES_CH[idx] = int(TES_charging.sum() > 0)
                TES_DC[idx] = int(TES_discharging.sum() > 0)
                
                # log output means
                gen[mean_idx]       = np.mean( pw.get_array('gen')[out_slice] )
                qdot_nuc[mean_idx]  = np.mean( pw.get_array('q_dot_nuc_inc')[out_slice] )
                defocus[mean_idx]   = np.mean( pw.get_array('defocus')[out_slice] )
                
                # log output std dev
                gen[stdv_idx]       = np.std( pw.get_array('gen')[out_slice] )
                qdot_nuc[stdv_idx]  = np.std( pw.get_array('q_dot_nuc_inc')[out_slice] )
                defocus[stdv_idx]   = np.std( pw.get_array('defocus')[out_slice] )
                
                del  op_mode_array, op_mode_str_profile
            
            # =========================================
            # Simulation failed! Log results
            # =========================================
            else:
                
                # check if last crash was due to Pyomo
                if dispatch:
                    pass
                
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
Storage['pyomo_bad_log']      = pyomo_bad_log
Storage['pyomo_bad_idx_log']  = pyomo_bad_idx_log

# locating output directory
output_dir = FileMethods.output_dir
filename = 'failureModes__{0}__2021_11__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_[{4},{5}]__PC_[{6},{7}].nuctes'.format(
                json, dispatch, sscH, pyoH, tshours.min(), tshours.max(), p_cycle.min(), p_cycle.max() )

NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'wb') as f:
    pickle.dump(Storage, f)

