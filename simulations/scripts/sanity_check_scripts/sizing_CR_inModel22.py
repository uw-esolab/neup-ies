#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 11:03:53 2021

@author: gabrielsoto
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import rc
import matplotlib.cm as cm
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)

from util.FileMethods import FileMethods
from scipy.optimize import curve_fit
import pint, pandas, os
u = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)

# known values
N   = 950*u.MW # full nuclear thermal power
eta = 0.41999  # power cycle efficiency

# guessed values
h  = 6*u.hr # equivalent full-load thermal storage hours
fN = np.linspace(0.5,0.9,100)   # fraction of full nuclear power going into power cycle @ design point
lN = np.linspace(0.5,0.8,100)    # fraction of thermal storage contributed to by nuclear energy
ff, ll = np.meshgrid(fN,lN)

# backing out some results
S  = h * ff * N / ll   # total amount of thermal energy storage available
Y  = S / h             # power cycle thermal rating
X  = eta * Y           # power cycle eletric output @ design point

# floating result
fR  = Y - ff*N          # full receiver thermal power


# =============================================================================
# Tower and Receiver and Heliostatfield Opt
# =============================================================================
filepath = os.path.join( FileMethods.data_dir, "Model2.2__tower_receiver_field_opt.csv" )
dataframe = pandas.read_csv(filepath) 

P_ref    = dataframe.loc[:,'P_ref (MWe)']
LandAcre = dataframe.loc[:,'land_area_base (acre)']
N_hel    = dataframe.loc[:,'N_hel']
h_rec    = dataframe.loc[:,'rec_height (m)']
D_rec    = dataframe.loc[:,'D_rec (m)']
h_tower  = dataframe.loc[:,'h_tower (m)']

P_ref_input = np.linspace( P_ref.min(), P_ref.max(), 1000 )


linear = lambda x,a,b: a*x+b

pNhel, pcov = curve_fit(linear, P_ref, N_hel)
N_hel_fit = lambda x: linear(x, *pNhel)

ph_rec, pcov = curve_fit(linear, P_ref, h_rec)
h_rec_fit = lambda x: linear(x, *ph_rec)

pd_rec, pcov = curve_fit(linear, P_ref, D_rec)
D_rec_fit = lambda x: linear(x, *pd_rec)


# SAM results
fig = plt.figure()
ax  = fig.gca()
axy = ax.twinx()
ax.plot( P_ref, N_hel, 'C2o-')
ax.plot( P_ref_input, N_hel_fit(P_ref_input), 'C2--')
axy.plot(P_ref, h_rec, 'C0o-', label='Height of Receiver')
axy.plot( P_ref_input, h_rec_fit(P_ref_input), 'C0--')
axy.plot(P_ref, D_rec, 'C1o-', label='Diameter of Receiver')
axy.plot( P_ref_input, D_rec_fit(P_ref_input), 'C1--')

ax.set_xlabel('PC Electric Output due to Solar (MWe) ', fontweight='bold')
ax.set_ylabel('# of Heliostats ', fontweight='bold')
axy.set_ylabel('Meters ', fontweight='bold')

axy.legend()

# =============================================================================
# Plots
# =============================================================================

def fmt_T(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} MWt" if plt.rcParams["text.usetex"] else f"{s} $MW_t$"

fig = plt.figure(figsize=(10,8))
ax = fig.gca()

# plotting 2D heat map for Qdot losses, with interpolation
im = ax.imshow(fR.m, origin='lower', interpolation='bicubic', extent=[fN.min(), fN.max(), lN.min(), lN.max()])

# creating colorbar for the 2D heatmap with label
cb = fig.colorbar(im, ax=ax)
cb.set_label(r'Full Receiver Thermal Output to PC @ Design Point ( MWt )', fontweight = 'bold')

#===================================
### First Contour: tank diameter

d_contour = ax.contour(fN, lN, X.m, origin='lower', colors='w')

# labels for each contour level
d_fmt = {}
for l in d_contour.levels:
    d_fmt[l] = r'PC : {0:.1f} $MW_e$'.format(l)
    
# Add contour labels
clb = ax.clabel(d_contour,fmt=d_fmt, fontsize=10)

#===================================
### First Contour: tank diameter

d_contour = ax.contour(fN, lN, fR.m, origin='lower', colors='C3')

# labels for each contour level
d_fmt = {}
for l in d_contour.levels:
    d_fmt[l] = r'Rec : {0:.1f} $MW_t$'.format(l)
    
# Add contour labels
clb = ax.clabel(d_contour,fmt=d_fmt, fontsize=10)

#===================================
### First Contour: tank diameter
# contour plot for TES cost
c_contour = ax.contour(fN, lN, N_hel_fit(fR.m*eta), colors='k') 

# labels for each contour level
c_fmt = {}
for l in c_contour.levels:
    c_fmt[l]= r'# of HS' 
    c_fmt[l]+=' : {0:.0f}'.format(l)
# Add contour labels
ax.clabel(c_contour,fmt=c_fmt, fontsize=10)


# setting x and y labels
ax.set_xlabel(r'$f_N$ - fraction of Nuclear thermal power to PC', fontweight = 'bold')
ax.set_ylabel(r'$\lambda_N$ - fraction of TS attributed to Nuclear', fontweight = 'bold')




def fmt_E(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} MWe" if plt.rcParams["text.usetex"] else f"{s} $MW_e$"


# Power Cycle Rating
fig = plt.figure(figsize=(10,8))
ax = fig.gca()
CS = ax.contour( fN,lN,fR.m*eta , origin='lower')
ax.clabel(CS, CS.levels, inline=True, fmt=fmt_E, fontsize=10)
ax.set_xlabel(r'$f_N$ - fraction of Nuclear thermal power to PC', fontweight='bold')
ax.set_ylabel(r'$\lambda_N$ - fraction of TS attributed to Nuclear', fontweight='bold')
ax.set_title('Power Cycle Electric Output for Receiver', fontweight='bold')