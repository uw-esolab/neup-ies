#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 00:30:11 2022

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
from util.FileMethods import FileMethods
import os, pint, time, copy, pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patheffects as PathEffects
from matplotlib.gridspec import GridSpec
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

jsons = ['model1_Hamilton_560_tariffx1',
         'model1_Hamilton_560_tariffx2',
         'model1_CAISO_Hamilton'
         ]

short_jsons = ['SAM Tariff Peaks x1',
         'SAM Tariff Peaks x2',
         'CAISO Tariffs'
         ]

# locating output directory
# output_dir = FileMethods.output_dir
output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper")

sensitivityStudy = "TESCost"  
start_name     = 'paramSweep' # testDefocus # failureModes #paramSweep
PySAM_name     = 'PySAM'  # PySAM  # ''
add_extra_Name = True
extra_name     = '2022_03'  # 2021_10  # 2021_11 # 2021_12
dispatch       = True # True # False
sscH           = 24   # 12 # 24
pyoH           = 48   # 24 # 48

# reference PPA value to compare for given tariff rate
compare = [6.537 , 6.58, 6.26]

# upper and lower limits for relative PPA values
plower  = [0.995, 0.875, 0.80]
pupper  = [1.04,  1.125, 1.05 ]

# optimal values
Pref_opt = [2, 1, 1]
TES_opt  = [0, 1, 1]

# ones array used to plot vertical lines later
vert = np.ones(500)

# generate name of file
def updated_json( json ):
    filename = '{0}_{1}__{2}__{3}__pyomo_{4:.0f}__horizon_{5:.0f}_{6:.0f}__TES_[{7},{8}]__PC_[{9},{10}]__sens_{11}.nuctes'.format(
                    start_name,
                    PySAM_name,
                    json,
                    extra_name if add_extra_Name else '',
                    dispatch, sscH, pyoH,
                    tshours.min(), tshours.max(),
                    p_cycle.min(), p_cycle.max(), sensitivityStudy )
    return filename

json_range = range(3)

fig = plt.figure(figsize=(14,8))
gs  = GridSpec(11,3, figure=fig)

for r in json_range:
    
    json_name = jsons[r] 
    
    if json_name is "model1_CAISO_Hamilton":
        # tshours   = np.array([    2,   4,   6 ])
        # p_cycle   = np.array([ 1100, 950, 800 ]) 
        tshours   = np.array([    3,   5,   7 ])
        p_cycle   = np.array([ 900, 750, 600 ]) 
        
    elif json_name is "model1_Hamilton_560_tariffx1":
        # tshours   = np.array([    0,   1,   2 ])
        # p_cycle   = np.array([  550, 500, 450 ]) 
        tshours   = np.array([    0,   1,   2 ])
        p_cycle   = np.array([  550, 500, 450 ]) 
        
    elif json_name is "model1_Hamilton_560_tariffx1_5" or "model1_Hamilton_560_tariffx2":
        # tshours   = np.array([    4,   6,   8 ])
        # p_cycle   = np.array([  950, 800, 650 ]) 
        tshours   = np.array([    3,   5,   7 ])
        p_cycle   = np.array([  850, 700, 550 ]) 
    
    # =============================================================================
    # extract data
    # =============================================================================
    
    filename = updated_json( json_name )
    
    # generate full path to file
    NTPath = os.path.join(output_dir, filename)
    
    # pickling
    with open(NTPath, 'rb') as f:
        Storage = pickle.load(f)
        
    # storage dictionary
    sensitivityStudy = Storage['sensitivityStudy']
    iterator0    = Storage['iterator0'] 
    iterator1    = Storage['iterator1']  
    iterator2    = Storage['iterator2']  
    ppa          = Storage['ppa'] 
    
    
    # setting figure title
    full_title = "PySAM - {0} tariffs {1} Pyomo - {2:.0f}/{3:.0f}hr horizons ".format(
                        json_name, "with" if dispatch else "without", sscH, pyoH)
    
    # =============================================================================
    # 
    # =============================================================================
    
    pcycle_plotrange = range(3)
    tscost_plotrange = range(5)
    colors = ['blueviolet', 'k', 'C5', 'red', 'orange']
    for i in pcycle_plotrange:
        
        array = copy.deepcopy( ppa[:,:,i] )
        array /= compare[r]
        slic = slice(4*r,4*r+3,1)
        
        if i == Pref_opt[r]:
            ax  = fig.add_subplot(gs[slic,i], facecolor='papayawhip')
        else:
            ax  = fig.add_subplot(gs[slic,i])
        if r == 2:
            ax.set_xlabel( 'tshours', fontsize=12, fontweight='bold' )
        
        ax.set_title( "P_cycle = {0} ".format(iterator2[i]), fontweight='bold' )
        
        if i == 0:
            ax.set_ylabel( "PPA Relative to Ref \n {0}".format(short_jsons[r]), fontsize=12,  fontweight='bold' )
        
        
        tol = 0.3
        t_margin = tshours[-1]*tol
        TES_opt_rnge = np.linspace( tshours[0]-t_margin, tshours[-1]+t_margin, len(vert)  )
        ax.plot(TES_opt_rnge, vert, '--k', label='Reference PPA')
        
        if i == Pref_opt[r]:
            TES_opt_array = tshours[TES_opt[r]] * vert
            Pref_opt_rnge = np.linspace( plower[r], pupper[r], len(vert)  )
            
            ax.plot( TES_opt_array, Pref_opt_rnge, '-.k', linewidth=2, label='Optimum Design at $28.4/kWh')
        
        for j,c in zip(tscost_plotrange, colors):
            
            ls = '-'
            
            if r == 0 and i == 2:
                ax.plot(iterator1, array[j,:],      color=c, linewidth=2, linestyle=ls,
                        label="TES cost = ${0}/kWh".format(iterator0[j]) )
                
                ax.legend(loc='best', fontsize=10, bbox_to_anchor=(1, 0.2),
                          fancybox=True, shadow=True)
            
            if c == 'k':
                ax.plot(iterator1, array[j,:], linestyle=ls, color=c, linewidth=2 )
            else:
                ax.plot(iterator1, array[j,:], linestyle=ls, color=c, linewidth=2 )
            
            ax.plot(iterator1, array[j,:], 'o', color='k', markersize=7 )
            ax.plot(iterator1, array[j,:], 'o', color=c,   markersize=5 )
            
            if (i,r) == (0,0):
                xtol    = 0.1
                xmargin = tshours[-1]*xtol
            ax.set_xlim( [tshours[0]-xmargin, tshours[-1]+xmargin ] )
            ax.set_ylim( [plower[r], pupper[r] ] )
            

plt.tight_layout()