#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 17:48:14 2022

@author: gabrielsoto
"""


import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy import interpolate
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

# =============================================================================
# json file names
# =============================================================================

jsons = ['model1_Hamilton_560_tariffx1',
         'model1_Hamilton_560_tariffx1_5',
         'model1_Hamilton_560_tariffx2',
         'model1_CAISO_Hamilton'
         ]

titles = ['Tariff Peaks x1',
          'Tariff Peaks x1.5',
          'Tariff Peaks x2',
          'CAISO']

nPlots = len(titles)


# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper")

# json name parameters
start_name     = 'paramSweep_varTurbineCost' # failureModes
PySAM_name     = 'PySAM' 
add_extra_Name = True
extra_name     = '2022_03'  #  2022_02
dispatch       = True 
sscH           = 24   
pyoH           = 48  
TES_min        = 0   
TES_max        = 8    # 7
PC_min         = 450  
PC_max         = 900  

# selecting coefficient array
coeffs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  

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
# for n,coeff in enumerate(coeff_list):
#     extr_str += "_{0}{1}".format( coeff.split('_')[1], str(coeffs[n]))

# generate name of file
def updated_json( json ):
    filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__{11}.nuctes'.format(
                    start_name,
                    PySAM_name,
                    json,
                    extra_name if add_extra_Name else '',
                    dispatch, sscH, pyoH,
                    TES_min, TES_max,
                    PC_min, PC_max, extr_str )
    return filename

# =============================================================================
# Loop through jsons
# =============================================================================

plotrange = range(nPlots)

timings_log = []
for i in plotrange:

    filename = updated_json( jsons[i] )
    
    # generate full path to file
    NTPath = os.path.join(output_dir, filename)
    
    # pickling
    with open(NTPath, 'rb') as f:
        Storage = pickle.load(f)
        
    # storage dictionary
    slicing = slice(0,15,1)
    
    timings  = Storage['time_elapsed'] 
    ppa      = Storage['ppa'][:,slicing] 
    
    # amount of comp time per full year sim
    N_sims = len( ppa.flatten() )
    time_per_dispatch_call = timings / N_sims 
    
    # amount of comp time per MILP+SAM calls (365 calls in 1 year)
    timings_log.append( time_per_dispatch_call.to('s').m / 365 ) 
    
timings_log = np.array( timings_log )
print( 'Average time elapsed per MILP+SAM call = {0} s'.format( timings_log.mean() ) )