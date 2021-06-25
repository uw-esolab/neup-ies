#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:05:11 2021

@author: gabrielsoto
"""


import modules.NuclearTES as NuclearTES
import os
import matplotlib.pyplot as plt
import numpy as np
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
from util.SSCHelperMethods import SSCHelperMethods

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

for i,th in enumerate(iterator1):
    for j,N in enumerate(iterator2):
        
        # print current position in loop
        print("tshours :         ", th)
        print("number of tanks : ", N)
        
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

        


# =============================================================================
# Plots
# =============================================================================

fig1 = plt.figure()
ax1 = fig1.gca()
im1 = ax1.imshow(dtank_array.T, origin='lower')
cb1 = plt.colorbar(im1)
cb1.set_label(r'Diameter (m)')
ax1.set_xlabel('tshours')
ax1.set_ylabel('# of Tanks')
ax1.set_xticks(range(len(iterator1.m)))
ax1.set_xticklabels(iterator1.m)
ax1.set_yticks(range(len(iterator2)))
ax1.set_yticklabels(iterator2)
plt.tight_layout()

fig2 = plt.figure()
ax2 = fig2.gca()
im2 = ax2.imshow(qdotLoss_array.T, origin='lower')
cb2 = plt.colorbar(im2)
cb2.set_label(r'$\dot{q}^{L}_{TES}$ (MW)')
ax2.set_xlabel('tshours')
ax2.set_ylabel('# of Tanks')
ax2.set_xticks(range(len(iterator1.m)))
ax2.set_xticklabels(iterator1.m)
ax2.set_yticks(range(len(iterator2)))
ax2.set_yticklabels(iterator2)
plt.tight_layout()

# wiPyomo24 = tuple( [0,range(len(iterator1)),range(len(iterator2))] )
# wiPyomo48 = tuple( [1,range(len(iterator1)),range(len(iterator2))] )
# woPyomo   = tuple( [2,range(len(iterator1)),range(len(iterator2))] )

# x_array = P_ref*iterator2

# array_list = [ ppa_array,
#                lcoe_nom_array ]

# colormarkers = ['C0', 'C1']

# label_list = [ 'PPA Price',
#                'LCOE Nom' ]

# lw = 2

# # Energy Outputs
# fig = plt.figure()
# ax = fig.gca()
# ax.plot( x_array, annual_energy_array[woPyomo] , 'C0',    linewidth = lw,  label='No Pyomo'  )
# ax.plot( x_array, annual_energy_array[wiPyomo48] , 'C0--',  linewidth = lw,  label='w/ Pyomo 48 hr'  )
# ax.plot( x_array, annual_energy_array[wiPyomo24] , 'C0:',  linewidth = lw,  label='w/ Pyomo 24 hr'  )

# ax.set_xlabel('P_ref', fontweight='bold')
# ax.set_ylabel('Annual Energy (TWh)', fontweight='bold')
# ax.legend(loc='best')
# ax.set_title('tshours = {0} hr'.format(iterator1[0]), fontweight='bold')

# # Financial Outputs
# fig = plt.figure()
# ax = fig.gca()

# for array, color, label in zip(array_list, colormarkers, label_list):
#     colorP1 = color + '--'
#     colorP2 = color + ':'
    
#     labelP1 = label + ' w/ Pyomo 48 Hr'
#     labelP2 = label + ' w/ Pyomo 24 Hr'
    
#     ax.plot( x_array, array[woPyomo]   , color,   linewidth = lw, label = label)
#     ax.plot( x_array, array[wiPyomo48] , colorP1, linewidth = lw, label = labelP1)
#     ax.plot( x_array, array[wiPyomo24] , colorP2, linewidth = lw, label = labelP2)

# ax.set_xlabel('P_ref', fontweight='bold')
# ax.set_ylabel('Price (¢/kWh)', fontweight='bold')
# ax.legend(loc='best')
# ax.set_title('tshours = {0} hr'.format(iterator1[0]), fontweight='bold')


