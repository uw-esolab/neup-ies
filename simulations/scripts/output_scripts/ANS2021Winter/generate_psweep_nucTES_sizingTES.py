#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:05:11 2021

@author: gabrielsoto
"""


import modules.NuclearTES as NuclearTES
import os, time
import matplotlib.pyplot as plt
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
from util.SSCHelperMethods import SSCHelperMethods
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

pid = os.getpid()
print("PID = ", pid)


# =============================================================================
#     Set Up
# =============================================================================

nuctes = NuclearTES.NuclearTES(is_dispatch=True)
u = nuctes.u

T_tes_hot  = (nuctes.SSC_dict['T_tank_hot_init']  * u.degC ).to('degK')
T_tes_cold = (nuctes.SSC_dict['T_tank_cold_init'] * u.degC ).to('degK')
tshours       = nuctes.SSC_dict['tshours'] * u.hr
q_nuc_thermal = nuctes.SSC_dict['q_dot_nuclear_des'] * u.MW
h_min   = nuctes.SSC_dict['h_tank_min'] * u.m
h_tank  = nuctes.SSC_dict['h_tank'] * u.m
n_tanks = nuctes.SSC_dict['tank_pairs']
u_tank  = nuctes.SSC_dict['u_tank'] * u.W / (u.m**2 * u.K)

# =============================================================================
#    Run Simulation
# =============================================================================

# setting up arrays to cycle through
tshours   = np.array([ 5., 7.5, 10., 12.5, 15.]) * u.hr
num_tanks = np.array([  1,   5,  10,   15, 20])

iterator1 = tshours
iterator2 = num_tanks

empty = np.zeros([  len(iterator1), len(iterator2) ])

dtank_array    = empty.copy()
qdotLoss_array = empty.copy() 
TEScost_array  = empty.copy() 

tic = time.perf_counter()
for i,th in enumerate(iterator1):
    for j,N in enumerate(iterator2):
        
        # print current position in loop
        print("tshours :         ", th)
        print("number of tanks : ", N)
        
# =============================================================================
#         run sim
# =============================================================================
        # defining directories
        nuctes = NuclearTES.NuclearTES(is_dispatch=True)
        
        # have to redo this step from the init
        nuctes.dispatch_wrap = nuctes.create_dispatch_wrapper( nuctes.PySAM_dict )
        
        # update SSC dictionary parameters
        nuctes.SSC_dict['tshours'] = th.m
        nuctes.SSC_dict['tank_pairs'] = N
        
        # run simulation
        nuctes.run_sim( run_loop=True )
        nt = nuctes.Plant
        so = nuctes.SO
        
# =============================================================================
#         estimates
# =============================================================================
        
        # defining Temperatures, Rho, Cp
        T_tes_ave = 0.5*(T_tes_hot + T_tes_cold)
        rho_ave   = SSCHelperMethods.get_rho_htf(u, T_tes_ave, 17)
        cp_ave    = SSCHelperMethods.get_cp_htf( u, T_tes_ave, 17)
        
        # Q dot to TES
        Q_tes_des = (q_nuc_thermal * th).to('MJ')
        
        # Volume of Tanks
        vol_one_temp_avail = Q_tes_des / (rho_ave * cp_ave * (T_tes_hot - T_tes_cold) )
        vol_one_temp_total = vol_one_temp_avail / (1.0 - h_min / h_tank)
        
        # Total Area of Tanks
        A_cs   = vol_one_temp_total / (h_tank * N) 
        
        # Diameter of each Tank
        d_tank = 2.0*(A_cs / np.pi)**(0.5) 
        d_tank = d_tank.to_compact()
        
        # Heat Transfer Losses of Tank
        UA_tank = u_tank*(A_cs + np.pi*d_tank*h_tank)*N
        q_dot_loss_des = UA_tank*(T_tes_ave - 15.0*u.degK)
        q_dot_loss_des = q_dot_loss_des.to_compact()
        
        print("--diameter:   ", d_tank)
        print("--q_dot_loss: ", q_dot_loss_des)

        # log outputs
        dtank_array[i,j]    = d_tank.to('m').m
        qdotLoss_array[i,j] = q_dot_loss_des.to('MW').m
        TEScost_array[i,j]  = nt.Outputs.csp_pt_cost_storage * 1e-6

        

toc = time.perf_counter()
print('Made it past the loop in ')
# =============================================================================
# Plots
# =============================================================================

# creating figure object
fig = plt.figure(constrained_layout=True)
ax = fig.gca()

#===================================
# plotting 2D heat map for Qdot losses, with interpolation
im = ax.imshow(qdotLoss_array.T, origin='lower', interpolation='bicubic')

# creating colorbar for the 2D heatmap with label
cb = fig.colorbar(im, ax=ax)
cb.set_label(r'$\dot{q}^{L}_{TES}$ (MW)', fontweight = 'bold')

# setting tick marks for x and y axes
ax.set_xticks(range(len(iterator1.m)))
ax.set_xticklabels(iterator1.m)
ax.set_yticks(range(len(iterator2)))
ax.set_yticklabels(iterator2)

# getting extent of heatmap to help with contours
xmin,xmax,ymin,ymax = im.get_extent()
x_series = np.linspace(xmin, xmax, len(iterator1.m))
y_series = np.linspace(ymin, ymax, len(iterator2))

#===================================
### First Contour: tank diameter
d_levels = [25, 45, 60, 100] # Define levels 
d_cmap = cm.binary

# contour plot for tank diameter
d_contour = plt.contour(x_series, y_series, dtank_array.T, d_levels, colors='w')

# labels for each contour level
d_fmt = {}
for l,s in zip(d_contour.levels, d_levels):
    d_fmt[l] = 'Diameter: ' + str(s) + ' m'
    
# Add contour labels
plt.clabel(d_contour,fmt=d_fmt)

#===================================
### Second Contour: TES cost
c_levels = [150, 200, 250, 300] # Define levels 
c_cmap = cm.binary

# contour plot for TES cost
c_contour = plt.contour(x_series, y_series, TEScost_array.T, c_levels, colors='k') 

# labels for each contour level
c_fmt = {}
for l,s in zip(c_contour.levels, c_levels):
    c_fmt[l] = '$' + str(s) + ' M'
    
# Add contour labels
plt.clabel(c_contour,fmt=c_fmt)


# setting x and y labels
ax.set_xlabel('tshours', fontweight = 'bold')
ax.set_ylabel('# of Tank Pairs', fontweight = 'bold')





# fig1 = plt.figure(constrained_layout=True)
# spec1 = gridspec.GridSpec(ncols=6, nrows=7, figure=fig1)
# ax1 = fig1.add_subplot(spec1[0:3, 0:3])
# ax2 = fig1.add_subplot(spec1[4:7, 0:3])
# ax3 = fig1.add_subplot(spec1[2:5, 3:6])


# im1 = ax1.imshow(dtank_array.T, origin='lower')
# cb1 = fig1.colorbar(im1, ax=ax1)
# cb1.set_label(r'Diameter (m)', fontweight = 'bold')
# ax1.set_xlabel('tshours', fontweight = 'bold')
# ax1.set_ylabel('# of Tanks', fontweight = 'bold')
# ax1.set_xticks(range(len(iterator1.m)))
# ax1.set_xticklabels(iterator1.m)
# ax1.set_yticks(range(len(iterator2)))
# ax1.set_yticklabels(iterator2)
# # plt.tight_layout()

# # fig2 = plt.figure()
# # ax2 = fig2.gca()
# im2 = ax2.imshow(qdotLoss_array.T, origin='lower')
# cb2 = fig1.colorbar(im2, ax=ax2)
# cb2.set_label(r'$\dot{q}^{L}_{TES}$ (MW)', fontweight = 'bold')
# ax2.set_xlabel('tshours', fontweight = 'bold')
# ax2.set_ylabel('# of Tanks', fontweight = 'bold')
# ax2.set_xticks(range(len(iterator1.m)))
# ax2.set_xticklabels(iterator1.m)
# ax2.set_yticks(range(len(iterator2)))
# ax2.set_yticklabels(iterator2)
# # plt.tight_layout()

# im3 = ax3.imshow(TEScost_array.T, origin='lower')
# cb3 = fig1.colorbar(im3, ax=ax3)
# cb3.set_label(r'TES Cost from SAM ($ million)', fontweight = 'bold')
# ax3.set_xlabel('tshours', fontweight = 'bold')
# ax3.set_ylabel('# of Tanks', fontweight = 'bold')
# ax3.set_xticks(range(len(iterator1.m)))
# ax3.set_xticklabels(iterator1.m)
# ax3.set_yticks(range(len(iterator2)))
# ax3.set_yticklabels(iterator2)
# plt.tight_layout()
