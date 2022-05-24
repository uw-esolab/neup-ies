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

# =============================================================================
# Mixing vs Non-mixing indirect charging
# =============================================================================

# empty lists to populate
mm = []
cc = []
p_mixwCSP = []
p_onlyLFR = []

# fractional mdot or qdot relative to actual oversized turbine output
ff = np.arange(0.4, 1.1, 0.1)

# a finer array of fractional qdot for plotting
qq = np.linspace(0,1.5,1000)

# looping over fractional heat input
for f in ff:
    
    # calculating slope and y-intercept of fractional power curve
    m,c = SSCHelperMethods.linearize_indirectTES_eff(0.77)
    mm.append(m)
    cc.append(c)
    
    # calculating fractional power relative to total output
    p_mixwCSP.append( m[0]*qq + c[0] )
    p_onlyLFR.append( m[1]*qq + c[1] )

# restructuring arrays
mm = np.array(mm)
cc = np.array(cc)
p_mixwCSP = np.array(p_mixwCSP)
p_onlyLFR = np.array(p_onlyLFR)


fig = plt.figure()
ax = fig.add_subplot(311)
ax.plot( ff, mm[:,0], 'C0', linewidth=3, label='Mix with CSP')
ax.plot( ff, mm[:,1], 'C1', linewidth=3, label='Just LFR')
ax.set_ylabel("Slope", fontweight='bold')
ax.legend(loc='best')

ax = fig.add_subplot(312)
ax.plot( ff, cc[:,0], 'C0', linewidth=3, label='Mix with CSP')
ax.plot( ff, cc[:,1], 'C1', linewidth=3, label='Just LFR')
ax.set_ylabel("Y-Intercept", fontweight='bold')
ax.legend(loc='best')

# =============================================================================
# plotting actual P vs Q
# =============================================================================

ind = 3
designP = 465 / ff[ind]
designQ = 950 / ff[ind]
P_mixwCSP = designP * p_mixwCSP[ind,:]
P_onlyLFR = designP * p_onlyLFR[ind,:]
Q_all = designQ * qq

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot( Q_all, P_mixwCSP, 'C0', linewidth=3, label='Mix with CSP')
ax.plot( Q_all, P_onlyLFR, 'C1', linewidth=3, label='Just LFR')
ax.set_ylabel("Power Curve", fontweight='bold')
ax.legend(loc='best')

# c_divide_f = np.array([c_/f_ for f_,c_ in zip(ff,cc)])
# nn = (465/950) * (mm + c_divide_f )

# # base linearization
# Pref    = 465 / ff

# # defining NE2 module
# dptes = DualPlantTES.DualPlantTES(is_dispatch=True)

# eta_p = []
# bb = []
# for P in Pref:
# # update SSC dictionary parameters
#     dptes.SSC_dict['P_ref'] = P
    
#     dptes.dispatch_wrap.set_design()
    
#     eta_p__, b__ = SSCHelperMethods.get_linearized_ud_params( dptes.ud_array, dptes.dispatch_wrap.q_pb_design, dptes.SSC_dict)
#     eta_p.append(eta_p__.m)
#     bb.append(b__.m)
    
# fig = plt.figure()
# ax = fig.add_subplot(311)
# ax.plot( ff, mm[:,0], linewidth=3, label='Mix with CSP')
# ax.plot( ff, mm[:,1], linewidth=3, label='Just LFR')
# ax.set_ylabel("Slope", fontweight='bold')
# ax.legend(loc='best')

# ax2 = fig.add_subplot(312)
# ax2.plot( ff, cc[:,0], linewidth=3)
# ax2.plot( ff, cc[:,1], linewidth=3)
# ax2.set_ylabel("Y-Intercept", fontweight='bold')

# ax3 = fig.add_subplot(313)
# ax3.plot( ff, nn[:,0], linewidth=3)
# ax3.plot( ff, nn[:,1], linewidth=3)
# ax3.set_ylabel("Linearized Eta", fontweight='bold')
# ax3.set_xlabel("Fraction of LFR / Oversized Turbine", fontweight='bold')

# fig = plt.figure()
# ax = fig.add_subplot(111)

# Q = np.linspace(100, 950*1.5, 1000)
# Pt = 465 / 0.8
# f_ind = np.where(ff - 0.8 > 0)[0][0]
# eta_HD = nn[f_ind,0]
# eta_LD = nn[f_ind,1]
# c_HD = cc[f_ind,0] * Pt
# c_LD = cc[f_ind,1] * Pt


# ax.plot(Q, Q * eta_LD  + c_LD, linewidth=3)
# ax.plot(Q, Q * eta_HD  + c_HD, linewidth=3)
