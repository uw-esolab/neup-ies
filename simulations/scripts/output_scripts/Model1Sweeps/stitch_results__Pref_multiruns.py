#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:06:06 2022

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
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
# Retrieving Data
# =============================================================================

# locating output directory
output_dir = FileMethods.output_dir

start_name     = 'failureModes' # testDefocus # failureModes
PySAM_name     = 'PySAM'  # PySAM  # ''
add_extra_Name = True
extra_name     = '2022_02'  # 2021_10  # 2021_11 # 2021_12
json_name      = 'model1_Hamilton_560_tariffx1_4'   # model1_CAISO_Hamilton  # model1_Hamilton_560_tariffx1_5  
dispatch       = True # True # False
sscH           = 24   # 12 # 24
pyoH           = 48   # 24 # 48
TES_min        = 0    # 0  # 2
TES_max        = 7   # 14
PC_min         = [ 1000,  800,  600,  450 ]  # 100 # 300 # 400 # 550
PC_max         = [ 1150,  950,  750,  550 ] # 500 # 850
# PC_min         = [900, 500, 200 ]  # 100 # 300 # 400 # 550
# PC_max         = [1200, 800, 400 ] # 500 # 850

# selecting coefficient array

coeffs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # wec tariffx2
# coeffs = [0.1, 1.0, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # tariffx1    // wecdsrx1
# coeffs = [0.3, 0.5, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # tariffx1_5  // wecdsrx1_5
# coeffs = [0.5, 0.7, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # tariffx2
# coeffs = [0.5, 0.7, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # CAISO
# coeffs = [1.0, 0.5, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # wec tariffx2
# coeffs = [0.5, 0.5, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # tariffx1.5 - BAD
# coeffs = [0.5, 0.2, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # tariffx2   - BAD


# key for coefficients
coeff_list = ['ec_p',     # proft term
              'ec_pcsu',  # PC cold SU term
              'ec_pchsu', # PC hot SU term 
              'ec_pcsd',  # PC SD term
              'ec_pcr',   # PC ramping term
              'ec_pcer',  # PC excess ramping term
              'ec_nsu',   # LFR cold SU termq
              'ec_nhsu',  # LFR hot SU term
              'ec_nsd',   # LFR SD term
              'ec_pcg',   # PC power generating term
              'ec_pcsb',  # PC SB term
              'ec_nt'     # LFR thermal power generating term
              ]

# coefficient strings
extr_str = ''
for n,coeff in enumerate(coeff_list):
    extr_str += "_{0}{1}".format( coeff.split('_')[1], str(coeffs[n]))


CombinedStorage = {}
for i, (min_, max_) in enumerate(zip(PC_min, PC_max)):

    # generate name of file
    filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__{11}.nuctes'.format(
                    start_name,
                    PySAM_name,
                    json_name,
                    extra_name if add_extra_Name else '',
                    dispatch, sscH, pyoH,
                    TES_min, TES_max,
                    min_, max_, extr_str )
    
    
    # generate full path to file
    NTPath = os.path.join(output_dir, filename)
    
    if os.path.exists( NTPath ):
        print("File found at {0}".format(NTPath) )
    else:
        print("File not found! Possible matches: ")
        file_list       = os.listdir( output_dir )
        relevant_sname  = '{0}_{1}__'.format( start_name, PySAM_name )
        relevant_filenames = np.array([f for f in file_list if relevant_sname in f])
        
        splits = np.array([ s.split(  relevant_sname )[1] for s in relevant_filenames ])
        splits = np.array([ s.split('__')[0] for s in splits])
        
        print_names = relevant_filenames[ splits == json_name ].tolist()
        for name in print_names:
            print( name )
        
    # =============================================================================
    # extract data
    # =============================================================================
    
    # pickling
    with open(NTPath, 'rb') as f:
        Storage = pickle.load(f)
    
    exceptions = ['tshours', 'sscH', 'pyoH', 'op_modes_list', 'turb_unit_cost', \
                  'tes_spec_cost', 'fin_yrs', 'fin_rate', 'json', 'dispatch' ]
    exceptions += coeff_list
    
    for key in Storage.keys():
        
        if i == 0:
            CombinedStorage[key] = Storage[key]
        else:
            if key not in exceptions:
                CombinedStorage[key] = np.hstack([ CombinedStorage[key], Storage[key]  ])



# locating output directory
filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__{11}.nuctes'.format(
                start_name,
                PySAM_name,
                json_name,
                extra_name if add_extra_Name else '',
                dispatch, sscH, pyoH,
                TES_min, TES_max,
                min(PC_min), max(PC_max), extr_str )
NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'wb') as f:
    pickle.dump(CombinedStorage, f)      
