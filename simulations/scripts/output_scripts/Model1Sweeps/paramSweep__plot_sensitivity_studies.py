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

short_jsons = ['m1_tariffx1',
         'm1_tariffx2',
         'm1_CAISO_'
         ]

# locating output directory
output_dir = FileMethods.output_dir
# output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper")

sensitivityStudy = "TESCost"  # HotTemp # TESCost
start_name     = 'paramSweep' # testDefocus # failureModes #paramSweep
PySAM_name     = 'PySAM'  # PySAM  # ''
add_extra_Name = True
extra_name     = '2022_02'  # 2021_10  # 2021_11 # 2021_12
json_name      = 'model1_Hamilton_560_tariffx2'   # model1_CAISO_Hamilton  # model1_Hamilton_560_tariffx1_5
dispatch       = True # True # False
sscH           = 24   # 12 # 24
pyoH           = 48   # 24 # 48

compare = [6.537 , 6.58, 6.26]

if sensitivityStudy == "TESCost":
    plower  = [0.98, 0.90, 0.80]
    pupper  = [1.05, 1.10, 1.10 ]
elif sensitivityStudy == "HotTemp":
    plower  = [0.98, 0.92, 0.86]
    pupper  = [1.07, 1.01, 0.96 ]

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

fig = plt.figure(figsize=(16,8))
gs  = GridSpec(11,3, figure=fig)

fig2 = plt.figure(figsize=(16,8))
gs2  = GridSpec(11,3, figure=fig2)
    

for r in json_range:
    
    json_name = jsons[r] 
    
    if json_name is "model1_CAISO_Hamilton":
        tshours   = np.array([    2,   4,   6 ])
        p_cycle   = np.array([ 1100, 950, 800 ]) 
        
    elif json_name is "model1_Hamilton_560_tariffx1":
        tshours   = np.array([    0,   1,   2 ])
        p_cycle   = np.array([  550, 500, 450 ]) 
        
    elif json_name is "model1_Hamilton_560_tariffx1_5" or "model1_Hamilton_560_tariffx2":
        tshours   = np.array([    4,   6,   8 ])
        p_cycle   = np.array([  950, 800, 650 ]) 
    
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
    
    plotrange = range(3)
    colors = ['C0', 'C1', 'C2']
    
    for i in plotrange:
        
        array = copy.deepcopy( ppa[:,i,:] )
        array /= compare[r]
        slic = slice(4*r,4*r+3,1)
        ax  = fig.add_subplot(gs[slic,i])
        if r == 2:
            ax.set_xlabel( sensitivityStudy )
        
        ax.set_title( "tshours = {0} ".format(iterator1[i]) )
        
        if i == 0:
            ax.set_ylabel( "PPA Relative to Ref \n {0}".format(short_jsons[r]) )
        
        for j,c in zip(plotrange, colors):
            
            if i == 0:
                ax.plot(iterator0, array[:,j],      color=c, linewidth=2, label="P_cycle = {0} ".format(iterator2[j]) )
                ax.legend(loc='best')
            else:
                ax.plot(iterator0, array[:,j],      color=c, linewidth=2 )
            ax.plot(iterator0, array[:,j], 'o', color=c )
            
            ax.set_ylim( [plower[r], pupper[r] ] )
            
    
    plotrange2 = range(5)
    colors = ['C0', 'C1', 'C2', 'C3', 'C4']
    for i in plotrange:
        
        array = copy.deepcopy( ppa[:,:,i] )
        array /= compare[r]
        
        ax2  = fig2.add_subplot(gs2[slic,i])
        if r == 2:
            ax2.set_xlabel( 'tshours' )
        
        ax2.set_title( "P_cycle = {0} ".format(iterator2[i]) )
        
        if i == 0:
            ax2.set_ylabel( "PPA Relative to Ref \n {0}".format(short_jsons[r]) )
        
        for j,c in zip(plotrange2, colors):
            
            if i == 0:
                ax2.plot(iterator1, array[j,:],      color=c, linewidth=2, label="{0} = {1} ".format(sensitivityStudy, iterator0[j]) )
                ax2.legend(loc='best')
            else:
                ax2.plot(iterator1, array[j,:],      color=c, linewidth=2 )
            ax2.plot(iterator1, array[j,:], 'o', color=c )
            
            ax2.set_ylim( [plower[r], pupper[r] ] )
            ax2.legend(loc='best')
                