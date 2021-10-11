#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 14:11:08 2021

@author: gabrielsoto
"""

from util.PySSCWrapper import PySSCWrapper
import os, cmath
from pylab import rc
import matplotlib.pyplot as plt
import pint
u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
import numpy as np
rc('axes', linewidth=2)
rc('font', weight='bold', size=12)

# print the PID of this script run
pid = os.getpid()
print("PID = ", pid)

# =============================================================================
# run SSC
# =============================================================================

# initialize the PySSC Wrapper
pw = PySSCWrapper(json_name='model2',is_debug=True)

# run SSC through PySSC
pw.run_sim()

# get output arrays from the SSC run
p_cycle       = pw.get_array('P_cycle') * u.MW
gen           = pw.get_array('gen') * u.kW
p_cool        = pw.get_array('P_cooling_tower_tot') * u.MW
q_dot_rec_in  = pw.get_array('q_dot_rec_inc') * u.MW
q_dot_nuc_in  = pw.get_array('q_dot_nuc_inc') * u.MW
q_startup     = pw.get_array('q_startup') * u.MW
m_dot_rec     = pw.get_array('m_dot_rec') * u.kg/u.hr
m_dot_nuc     = pw.get_array('m_dot_nuc') * u.kg/u.s
T_pc_in       = pw.get_array('T_pc_in') * u.degC
T_pc_out      = pw.get_array('T_pc_out') * u.degC
e_ch_tes      = pw.get_array('e_ch_tes')
defocus       = pw.get_array('defocus')
t_plot        = (pw.get_array('time_hr') / 24 ) *u.d
op_mode_1     = pw.get_array('op_mode_1')

q_ch_tes  = pw.get_array('q_ch_tes')
q_dc_tes  = pw.get_array('q_dc_tes')

T_nuc_in       = pw.get_array('T_nuc_in') * u.degC
T_rec_in       = pw.get_array('T_rec_in') * u.degC
T_nuc_out      = pw.get_array('T_nuc_out') * u.degC
T_rec_out      = pw.get_array('T_rec_out') * u.degC

T_tes_cold_in   = pw.get_array('T_tes_cold_in') * u.degC
T_tes_cold_out  = pw.get_array('T_tes_cold') * u.degC
T_tes_hot       = pw.get_array('T_tes_hot') * u.degC

# =============================================================================
# Temperature differences
# =============================================================================

# plotting temp out differences
fig = plt.figure()
ax = fig.gca()

label_CSP_off = True
label_CSP_on  = True
label_CSP_su  = True

log_OFF_times = []
log_ON_times = []

for t in range(len(T_rec_out)):
    T_out_diff = T_nuc_out-T_rec_out
    
    if q_dot_rec_in[t] == 0:
        log_OFF_times.append(t)
        if q_startup[t] > 0:
            if label_CSP_su:
                ax.plot(t_plot[t], T_out_diff[t], 'C3.', linewidth=2, label='LFR ON, CSP SU')
                label_CSP_su = False
            else:
                ax.plot(t_plot[t], T_out_diff[t], 'C3.', linewidth=2)
        else:
            if label_CSP_off:
                ax.plot(t_plot[t], T_out_diff[t], 'C1.', linewidth=2, label='LFR ON, CSP OFF')
                label_CSP_off = False
            else:
                ax.plot(t_plot[t], T_out_diff[t], 'C1.', linewidth=2)
    else:
        log_ON_times.append(t)
        if label_CSP_on:
            ax.plot(t_plot[t], T_out_diff[t], 'C0.', linewidth=2, label='LFR ON, CSP ON')
            label_CSP_on = False
        else:
            ax.plot(t_plot[t], T_out_diff[t], 'C0.', linewidth=2)
        
ax.set_xlabel('Time (d)', fontweight='bold')
ax.set_ylabel('Diff T_salt_out b/e Nuc and CSP \n(deg C)', fontweight='bold')
ax.legend(loc='best')


# =============================================================================
# enthalpies
# =============================================================================

def c_p( T_C ):
    """ c_p in units of kJ/kgK
    """
    T_K = T_C.to('degK').m
    T_C = T_C.m
    c_p = -1E-10*T_K*T_K*T_K + 2E-07*T_K*T_K + 5E-06*T_K + 1.4387
    c_p = c_p * u.kJ / u.kg / u.K
    
    # c_p = 1443. + 0.172 * (T_C)
    # c_p = c_p * u.J / u.kg / u.degC
    return c_p.to('kJ/kg/K')

def Temp( H ):
    """ H in units of J/kg
    """
    H = H.to('J/kg').m
    T_K = -0.0000000000262*H*H + 0.0006923 * H + 0.03058
    T_K = T_K*u.degK
    return T_K.to('degC')

def enth( T_C ):
    """ H in units of J/kg
    """
    T_C = T_C.to('degC').m
    H = 1443.*T_C + 0.086*T_C*T_C
    H = H * u.J/u.kg
    return H

# =============================================================================

Tc_range = np.linspace(290,600,1000)*u.degC
Cp_range = c_p(Tc_range)

diffC = Cp_range[-1] - Cp_range[0]
diffC /= Cp_range[0]
diffC = diffC.m

fig = plt.figure()
ax = fig.gca()
ax.plot(Tc_range, Cp_range, '.', label=r"Pct Diff = {:.2%}" .format(diffC) )
ax.set_xlabel('Temperature (deg C)', fontweight='bold')
ax.set_ylabel('Specific Heat of NaNO3-KNO3 \n (kJ/(kg K))', fontweight='bold')
ax.legend(loc='best')

# =============================================================================

a1, a2, a3, a4 = [-1E-10, 2E-07, 5E-06, 1.4387]
a1 *= u.kJ/u.kg/u.K**4
a2 *= u.kJ/u.kg/u.K**3
a3 *= u.kJ/u.kg/u.K**2
a4 *= u.kJ/u.kg/u.K

def intg_enthalpy(T):
    """ Calculate enthalpy from integrated Cp empirical formula 
    """
    T = T.to('degK')
    h = a1/4*T**4 + a2/3*T**3 + a3/2*T**2 + a4*T
    return h.to('kJ/kg')

# enthalpies from each stream node
h_in      = intg_enthalpy(T_nuc_in)
h_nuc_out = intg_enthalpy(T_nuc_out)
h_rec_out = intg_enthalpy(T_rec_out)
h_rec_out[log_OFF_times] *= 0

# change in enthalpy across each plant
dh_nuc = h_nuc_out - h_in
dh_rec = h_rec_out - h_in

# mass flow of mix
m_dot_nuc = m_dot_nuc.to('kg/s')
m_dot_rec = m_dot_rec.to('kg/s')
m_dot_mix = m_dot_nuc + m_dot_rec

# change in enthalpy across mixture -> through mass-flow weighted average
dh_mix = (dh_nuc*m_dot_nuc + dh_rec*m_dot_rec) / m_dot_mix


X = np.zeros([len(t_plot),4,2]) * u.degK
for i,t in enumerate(t_plot):
    
    A = a1/4
    B = a2/3
    C = a3/2
    D = a4
    E = -(h_in[i] + dh_mix[i])
    
    a = (B/A).to('degK').m
    b = (C/A).to('degK^2').m
    c = (D/A).to('degK^3').m
    d = (E/A).to('degK^4').m
    
    M = np.array([[0,0,0,-d],
                  [1,0,0,-c],
                  [0,1,0,-b],
                  [0,0,1,-a]])
    
    U,V = np.linalg.eig(M)
    U = U*u.degK
    X[i,0,0] = U[0].real
    X[i,0,1] = U[0].imag
    
    X[i,1,0] = U[1].real
    X[i,1,1] = U[1].imag
    
    X[i,2,0] = U[2].real
    X[i,2,1] = U[2].imag
    
    X[i,3,0] = U[3].real
    X[i,3,1] = U[3].imag

T_mix = X[:,2,0].to('degC')


# =============================================================================
# now need to solve a quartic formula for T of the mixture
# log_beta = []
# X = np.zeros([len(t_plot),4,2]) * u.degK
# for i,t in enumerate(t_plot):
#     A = a1/4
#     B = a2/3
#     C = a3/2
#     D = a4
#     E = -(h_in[i] + dh_mix[i])
    
#     alpha = -(3/8.) * (B**2/A**2) + C/A
#     beta  = B**3/(8*A**3) - B*C/(2*A**2) + D/A
#     gamma = -(3/256)*(B**4/A**4) + C*B**2/(16*A**3) - B*D/(4*A**2) + E/A
    
#     log_beta.append(beta==0)
    
#     P = -alpha**2/12 - gamma
#     Q = -alpha**3/108 + alpha*gamma/3 - beta**2/8
#     QParg = Q**2/4 + P**3/27
#     QPu   = QParg.u
#     R = -Q/2 + np.sqrt(QParg + 0j)
#     U = R**(1/3.)
    
#     y = -5/6*alpha + 0j
#     y +=  (-(Q.m)**(1/3) if U == 0 else (U - P/(3*U)).m )*y.u
    
#     W = (alpha + 2*y)**(1/2)
    
#     s_sign = [ 1,  1, -1, -1]
#     t_sign = [ 1 ,-1,  1, -1]
    
#     count = 0
#     for s,t in zip(s_sign,t_sign):
#         Xn = -B/(4*A) + (s*W + t*(-(3*alpha + 2*y + s*2*beta/W))**(1/2))/2
        
#         X[i,count,0] = Xn.real
#         X[i,count,1] = Xn.imag
        
#         count += 1

# T_mix = X[:,1,0].to('degC')

fig = plt.figure()
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
# ax = fig.gca()
ax1.plot(t_plot, T_mix, 'C0.', label='Mix')
ax1.plot(t_plot, T_nuc_out, 'C1.', label='LFR')
ax1.plot(t_plot, T_rec_out, 'C2.', label='CSP')
ax1.set_xlabel('Time (d)', fontweight='bold')
ax1.set_ylabel('Outlet Temperatures \n(deg C)', fontweight='bold')
ax1.legend(loc='best')

ax2.plot(t_plot, T_mix, '.', label='mix')
ax2.plot(t_plot, T_nuc_out, '.', label='lfr')
ax2.plot(t_plot, T_rec_out, '.', label='csp')
ax2.set_xlabel('Time (d)', fontweight='bold')

# =============================================================================
# now just getting the mass-flow weighted outlet temperatures
    
# change in enthalpy across mixture -> through mass-flow weighted average
T_avg_mix = (T_nuc_out*m_dot_nuc + T_rec_out*m_dot_rec) / m_dot_mix
T_avg_mix = T_avg_mix.to('degC')
# T_avg_mix[log_OFF_times] *= 0

# fig = plt.figure()
# ax = fig.gca()
# ax.plot(t_plot, T_mix.to('degC'), '.',  label='enthalpy averaged')
# ax.plot(t_plot, T_avg_mix.to('degC'), '.',  label='mass flow averaged')
# ax.set_xlabel('Time (d)', fontweight='bold')
# ax.set_ylabel('Outlet Temperatures Diff \n(deg C)', fontweight='bold')
# ax.legend(loc='best')

diffT = T_mix.to('degC')-T_avg_mix.to('degC')
diffT = diffT.to('mK')

fig = plt.figure()
ax = fig.gca()
ax.plot(t_plot, diffT, '.',  label='diff b/e int and avg')
# ax.plot(t_plot[log_ON_times], (T_mix-T_avg_mix).to('degC')[log_ON_times], '.',  label='CSP ON')
ax.set_xlabel('Time (d)', fontweight='bold')
ax.set_ylabel('Outlet Temperatures Diff \n(deg mK)', fontweight='bold')
ax.set_title('Difference between Enthalpy-Averaged and MassFlow-Averaged T_out_mix', fontweight='bold')
# ax.legend(loc='best')