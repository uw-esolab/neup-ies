#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 17:45:28 2022

@author: gabrielsoto
"""

from util.SSCHelperMethods import SSCHelperMethods
import numpy as np
import matplotlib.pyplot as plt
import modules.DualPlantTES as DualPlantTES
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)

# Mixing vs Non-mixing indirect charging
mm = []
cc = []

ff = np.arange(0.4, 1.0, 0.05)

for f in ff:
    m,c = SSCHelperMethods.linearize_indirectTES_eff(f)
    mm.append(m)
    cc.append(c)

mm = np.array(mm)
cc = np.array(cc)

c_divide_f = np.array([c_/f_ for f_,c_ in zip(ff,cc)])
nn = (465/950) * (mm+ c_divide_f )

# base linearization
Pref    = 465 / ff

# defining NE2 module
dptes = DualPlantTES.DualPlantTES(is_dispatch=True)

eta_p = []
bb = []
for P in Pref:
# update SSC dictionary parameters
    dptes.SSC_dict['P_ref'] = P
    
    dptes.dispatch_wrap.set_design()
    
    eta_p__, b__ = SSCHelperMethods.get_linearized_ud_params( dptes.ud_array, dptes.dispatch_wrap.q_pb_design, dptes.SSC_dict)
    eta_p.append(eta_p__.m)
    bb.append(b__.m)
    
fig = plt.figure()
ax = fig.add_subplot(311)
ax.plot( ff, mm[:,0], linewidth=3, label='Mixing')
ax.plot( ff, mm[:,1], linewidth=3, label='No Mixing')
ax.set_ylabel("Slope", fontweight='bold')
ax.legend(loc='best')

ax2 = fig.add_subplot(312)
ax2.plot( ff, cc[:,0], linewidth=3)
ax2.plot( ff, cc[:,1], linewidth=3)
ax2.set_ylabel("Y-Intercept", fontweight='bold')

ax3 = fig.add_subplot(313)
ax3.plot( ff, nn[:,0], linewidth=3)
ax3.plot( ff, nn[:,1], linewidth=3)
ax3.set_ylabel("Linearized Eta", fontweight='bold')
ax3.set_xlabel("Fraction of LFR / Oversized Turbine", fontweight='bold')

# ax4 = fig.add_subplot(414)
# ax4.plot( ff, eta_p, linewidth=3)
# ax4.set_ylabel("Linearized Eta (Original)", fontweight='bold')
# ax4.set_xlabel("Fraction of LFR / Oversized Turbine", fontweight='bold')