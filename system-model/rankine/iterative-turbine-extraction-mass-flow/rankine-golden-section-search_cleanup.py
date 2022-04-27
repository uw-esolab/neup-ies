from state import *
from components import *

import numpy as N

import time
from math import log10, ceil

# --------------------------
def log_header(N,names):
    hh = []
    for name in names:
        h = []
        for i in range(1,N):
            h.append("{:s}{:d}".format(name,i))
        hh.append(h)
    return ','.join([','.join(h) for h in hh]) + '\n'

def log_state(*args):
    # T,P,s,h,x,m
    return ','.join([','.join([str(v) for v in arr[1:]]) for arr in args]) + '\n'

step = lambda xold,xnew: xold + (xnew-xold)*0.66
# --------------------------

def update_mf(m_dot_new, m_dot_old):
    return m_dot_old + (m_dot_new - m_dot_old)*0.15

tstart = time.time()

DT = 5.55 #[K]    #Target Approach Temperature 
 
N_hxrs = 20 #[-]    #Number of heat exchangers

ns = 32  #number of states used + 1
 
fout = open('iter-log.csv','w')
fout.write(log_header(ns, "T P h x m".split(' ')))

# initialize arrays
m = N.zeros(ns)
T_H_A = N.zeros(3*N_hxrs+2)
T_C_A = N.zeros(3*N_hxrs+2)
T_H_B = N.zeros(3*N_hxrs+2)
T_C_B = N.zeros(3*N_hxrs+2)
T_H_C = N.zeros(3*N_hxrs+2)
T_C_C = N.zeros(3*N_hxrs+2)
T_H_D = N.zeros(3*N_hxrs+2)
T_C_D = N.zeros(3*N_hxrs+2)
T_H_E = N.zeros(3*N_hxrs+2)
T_C_E = N.zeros(3*N_hxrs+2)
T_H_F = N.zeros(3*N_hxrs+2)
T_C_F = N.zeros(3*N_hxrs+2)
T_H_G = N.zeros(3*N_hxrs+2)
T_C_G = N.zeros(3*N_hxrs+2)
h = N.zeros(ns)
P = N.zeros(ns)
T = N.zeros(ns)
m = N.zeros(ns)
m_t = N.zeros(ns)
x = N.zeros(ns)

loc = N.zeros(N_hxrs*3+2)

for i in range(1,(3*N_hxrs+1)+1):
    loc[i] = (i-1)/(3*N_hxrs)
 
 
 
PR = 2.657 #[-]  #even pressure accross turbines assumption
 
#Mass Flow Rates ==========
 
# "Turbines"
# {m[3] = 55{50} [kg/s]    #HPT outlet}
# {m[6] = 48{43.9} [kg/s]    #IPT1 outlet}
# {m[7] = 35{19.4} [kg/s]    #IPT2 outlet}
# m[8] = 28.7 #[kg/s]    #IPT3 outlet
# {m[10] = 22 [kg/s]    #LPT1 outlet}
# {m[11] = 18{12.99} [kg/s]    #LPT2 outlet}
# {m[12] = 15 {17.29} [kg/s]    #LPT3 outlet}
# {m[13] = 13.57 [kg/s]}    #LPT4 outlet
 
m[1] = 477.2 #[kg/s]
 
# ----------- initial calculations based on guesses
f_hp3 = 0.00214
f_ip6 = 0.1254
f_ip7 = 0.1026
f_ip8 = 0.0768
f_lp10 = 0.01797
f_lp11 = 0.05112
f_lp12 = 0.04444
f_lp13 = 0.03759
m[23] = m[1]
h[20] = 436116
h[23] = 881554
h[24] = 1146000
h[27] = 663252
h[28] = 420997
h[29] = 316028
h[30] = 226272
# ---------------*******

# Initial guess for T[1]
T[1] = 611.4 #[C]
 
it = 0
it_max = 50
err_iter = 999.
tol = 0.0001

while((err_iter > tol and it<it_max) or it<3 ):

    print("Iter {:d}: ".format(it), end='')

    m[2] = m[1]    #Boiler
    m[3] = m[2]*f_hp3    #HPT outlet fraction
    m[4] = m[3]    #Across CSP 
    m[5] = m[2] - m[3]    #HPT to IPT
    m[6] = m[5]*f_ip6    #IPT1 outlet fraction
    m[7] = (m[5] - m[6])*f_ip7    #IPT2 outlet fraction
    m[8] = (m[5] - m[6] - m[7])*f_ip8    #IPT3 outlet fraction
    m[9] = m[5] - m[6] - m[7] - m[8]    #IPT to LPT
    m[10] = m[9]*f_lp10    #LPT1 outlet fraction
    m[11] = (m[9]-m[10])*f_lp11    #LPT2 outlet fraction
    m[12] = (m[9] - m[10] - m[11])*f_lp12    #LPT3 outlet fraction
    m[13] = (m[9] - m[10] - m[11] - m[12])*f_lp13    #LPT4 outlet fraction
    m[14] = m[9] - m[10] - m[11] - m[12] - m[13]    #LPT to condenser
    # -----------------------------

    P[2] = 30e6     #Boiler outlet pressure
    T[2] = 632+273.15    #Boiler outlet temperature
    h[2] = enthalpy( T = T[2], P = P[2])    #Boiler outlet enthalpy

    #LFR Boiler -----------------------------------------
    delta_P_boiler = 5*286107  #[Pa]    #Assumption of pressure drop accross boiler (Jacob Wenner)
    P[1] = P[2] + delta_P_boiler    #Boiler inlet pressure
    
    # {T[1] = 340+273.15    #Boiler inlet temperature}
    h[1] = enthalpy( T = T[1], P = P[1])    #Boiler inlet enthalpy
    
    Q_dot_LFR = 950e6    #Heat transfer from LFR
    
    h[2] = (m[1]*h[1] + Q_dot_LFR)/m[2]    #Energy balance on Boiler
    

    # constants 
    eta_HPT = 0.90 #[-]     #HPT isentropic efficiency
    eta_IPT1 = 0.90 #[-]    #IPT1 isentropic efficiency
    eta_IPT2 = 0.90 #[-]     #IPT2 isentropic efficiency
    eta_IPT3 = 0.90 #[-]    #IPT3 isentropic efficiency
    eta_LPT1 = 0.90 #[-]     #LPT1 isentropic efficiency
    eta_LPT2 = 0.90 #[-]    #LPT2 isentropic efficiency
    eta_LPT3 = 0.90 #[-]    #LPT3 isentropic efficiency
    eta_LPT4 = 0.90 #[-]    #LPT4 eisentropic efficiency
    eta_LPT5 = 0.90 #[-]     #LPT5 eisentropic efficiency
    DT_int_G = 6.5 #[K]
    DT_G = 5.55 #[K]
    DT_int_F = 3 #[K]
    DT_F = 5.55 #[K]
    DT_int_E = 3 #[K]
    DT_E = 5.55 #[K]
    eta_HP = 0.85 #[-]    #HP pump efficiency
    DT_int_D = 6.5 #[K]
    DT_D = 5.55 #[K]
    DT_int_C = 3 #[K]
    DT_C = 5.55 #[K]
    DT_int_B = 3#[K]
    DT_B = 5.55 #[K]
    DT_int_A = 3 #[K]
    DT_A = 5.55 #[K]
    eta_LP = 0.85 #[-]    #LP pump efficiency
    T_r = 25+273.15    #condenser Reservoir Temperature
    
    P[16] = 17.24*1e5    #LP pump pressure outlet
    P[14] = 5000 #[Pa]    #LPT5 outlet pressure
    x[14] = 0.9 #[-]    #LPT5 outlet quality
    P[13] =P[14]*PR     #LPT4 outlet pressure
    P[12] = P[13]*PR    #LPT3 outlet pressure
    P[11] = P[12]*PR   #LPT2 outlet pressure
    P[10] =P[11]*PR   #LPT1 outlet pressure
    P[8] = P[10]*PR   #IPT3 outlet pressure
    P[7] = P[8]*PR    #IPT2 outlet pressure
    P[6] = P[7]*PR    #IPT1 outlet pressure
    P[3] = P[6]*PR    #HPT outlet pressure

    #HP Turbine -----------------------------------------
    
    # "HPT"
    h[3], W_HPT = turbine( m[2], h[2], P[2], P[3], eta_HPT)
    T[3] = temperature( P = P[3], H = h[3])    #HPT outlet temperature
    
    #IP Turbine -------------------------------------------
    
    # "IPT1"
    T[5] = T[3]    #IPT1 inlet temperature
    P[5] = P[3]    #IPT1 inlet pressure
    h[5] = h[3]    #IPT1 inlet enthalpy
    
    h[6],W_IPT1 = turbine( m[5], h[5], P[5], P[6], eta_IPT1)    
    # T[6] = temperature( H = h[6], P = P[6])    #IPT1 outlet temperature
    
    # "IPT2"
    m_t[6] = m[5]-m[6]    #IPT2 mass flow rate
    h[7], W_IPT2 = turbine( m_t[6], h[6],P[6],P[7],eta_IPT2)    
    # T[7] = temperature( P = P[7], H = h[7])    #IPT2 outlet temperature
    
    # "IPT3"    
    m_t[7] = m_t[6]-m[7]    #IPT3 mass flow rate
    h[8],W_IPT3 = turbine( m_t[7], h[7], P[7], P[8], eta_IPT3)    
    # T[8] = temperature( P = P[8], H = h[8])    #IPT3 outlet temperature
    
    
    #LP Turbine ------------------------------------------
    
    # "LPT1"
    # T[9] = T[8]     #LPT1 inlet temperature
    P[9] = P[8]    #LPT1 inlet pressure
    h[9] = h[8]    #LPT1 inlet enthalpy
    
    h[10], W_LPT1 = turbine( m[9], h[9],P[9],P[10],eta_LPT1)
    
    # T[10] = temperature( P = P[10], H = h[10])    #LPT1 outlet temperature
    
    # "LPT2"
    m_t[10] = m[9]-m[10]    #LPT2 mass flow rate
    
    h[11],W_LPT2 = turbine( m_t[10], h[10], P[10], P[11], eta_LPT2)
    
    # T[11] = temperature( P = P[11], H = h[11])     #LPT2 outlet temperature
    
    # "LPT3"
    m_t[11] = m_t[10]-m[11]    #LPT3 mass flow rate
    h[12],W_LPT3 = turbine( m_t[11], h[11], P[11], P[12], eta_LPT3)
    # T[12] = temperature( P = P[12], H = h[12])     #LPT3 outlet temperature
    
    # "LPT4"
    m_t[12] = m_t[11]-m[12]    #LPT4 mass flow rate
    h[13],W_LPT4 = turbine( m_t[12], h[12], P[12], P[13], eta_LPT4)
    # T[13] = temperature( P = P[13], H = h[13])    #LPT4 outlet temperature
    
    # "LPT5"
    m_t[13] = m_t[12] - m[13]     #LPT5 mass flow rate
    h[14],W_LPT5 = turbine( m_t[13], h[13], P[13], P[14], eta_LPT5)
    # T[14] = temperature( P = P[14], H = h[14])    #LPT5 outlet temperature
    
    #CSP -----------------------------------------------------
    
    # "temp equations for csp"
    T[4] = T[3]
    P[4] = P[3]
    h[4] = h[3]
    
    #HX G Closed FWH -------------------------------
    P[24] = P[1] 
    P[25] = P[4]
    m[24] = m[23]
    
    print("..HXG", end='')
    m4_req, h_test_1_, h[25], T_H_G[1:3*N_hxrs+1+1], T_C_G[1:3*N_hxrs+1+1] = fwh(0, 0, 0, m[24], h[24], P[24], h[4], P[4], DT_int_G, DT_G, N_hxrs)
    m[4] = update_mf(m4_req, m[4])
    m[3] = m[4]    #Across CSP 
    f_hp3 = m[3]/m[2]    #HPT outlet fraction
    m[5] = m[2] - m[3]    #HPT to IPT
    m[25] = m[4] 
    # T_test[1] = temperature( H = h_test[1], P = P[1])
    
    m[28] = m[10]
    m[29] = m[11]+m[28]
    m[30] = m[12]+m[29]
    m[31] = m[13]+m[30]
    m[15] = m[14]+m[31]    #Combine flow mass balance
    m[16] = m[15]
    m[17] = m[16]
    m[18] = m[17]
    m[19] = m[18]
    m[20] = m[19]
    #collected mass flow equations from elsewhere
    m[26] = m[25]+m[6]
    m[27] = m[26]+m[7]
    m[21] = m[20] + m[8] + m[27]
    m[22] = m[21]
    m[23] = m[22]
    #..

    P[17] = P[16]
    P[30] = P[12]
    P[18] = P[17]
    P[29] = P[11]
    P[19] = P[18]
    P[28] = P[10]
    P[20] = P[19]
    P[23] = P[24]
    P[26] = P[6] 
    P[22] = P[23]
    P[27] = P[7]

    #HX Open FWH --------------------------------------
    print("..OpenFWH", end='')
     
    # T[27] = temperature( P = P[27], H = h[27])
    m21_req, h[21], P[21] = open_fwh( m[20], h[20], P[20], m[8], h[8], P[8], m[27], h[27], P[27] )
    m[21] = update_mf(m21_req, m[21])

    m[22] = m[21]
    m[23] = m[22]
    # T[21] = temperature( H = h[21], P = P[21])    #Open FWH outlet temperature
    
    #High Pressure Pump -----------------------------
    h[22], W_dot_HP = pump( m[21], h[21], P[21], P[22], eta_HP)     
    
    
    #Condenser --------------------------------------------
    h[15], P[15], Q_dot_cond, epsilon_cond = condenser( m[15], h[14], P[14], T_r, DT)
    # T[15] = temperature( H = h[15], P = P[15])    #Condenser temperature outlet

    #Low Pressure Pump -----------------------------
    h[16], W_dot_LP = pump( m[15], h[15], P[15], P[16], eta_LP)
    # T[16] = temperature( P = P[16], H =h[16])

    #HX A Closed FWH ---------------------------------
    print("..HXA", end='')
    # T[30] = temperature( P = P[30], H = h[30])
    m13_req, h[17], h[31], T_H_A[1:3*N_hxrs+1+1], T_C_A[1:3*N_hxrs+1+1] = fwh( m[30], h[30], P[30], m[16], h[16], P[16], h[13], P[13], DT_int_A, DT_A, N_hxrs)
    m[13] = update_mf(m13_req, m[13])

    # T[17] = temperature( P = P[17], H = h[17])
    P[31] = P[13]

    #HX B Closed FWH ----------------------------------
    print("..HXB", end='')
    # T[29] = temperature( P = P[29], H = h[29])
    m12_req, h[18], h[30], T_H_B[1:3*N_hxrs+1+1], T_C_B[1:3*N_hxrs+1+1] = fwh( m[29], h[29], P[29], m[17], h[17], P[17], h[12], P[12], DT_int_B, DT_B, N_hxrs)
    m[12] = update_mf(m12_req, m[12])

    # T[18] = temperature( P = P[18], H = h[18])

    #HX C Closed FWH ---------------------------------
    print("..HXC", end='')
    # T[28] = temperature( P = P[28], H = h[28])
    m11_req, h[19], h[29], T_H_C[1:3*N_hxrs+1+1], T_C_C[1:3*N_hxrs+1+1] = fwh( m[28], h[28], P[28], m[18], h[18], P[18], h[11], P[11], DT_int_C, DT_C, N_hxrs)
    m[11] = update_mf(m11_req, m[11])
    # T[19] = temperature( P = P[19], H = h[19])

    #HX D Closed FWH ---------------------------------
    print("..HXD", end='')
    m10_req, h[20], h[28], T_H_D[1:3*N_hxrs+1+1], T_C_D[1:3*N_hxrs+1+1] = fwh( 0, 0, 0, m[19], h[19], P[19], h[10], P[10], DT_int_D, DT_D, N_hxrs)
    m[10] = update_mf(m10_req, m[10])
    # T[20] = temperature( P = P[20], H = h[20])


    #HX F Closed FWH --------------------------------
    print("..HXF", end='')
    # T[25] = temperature( P = P[25], H = h[25])
    m6_req, h[24], h[26], T_H_F[1:3*N_hxrs+1+1], T_C_F[1:3*N_hxrs+1+1] = fwh( m[25], h[25], P[25], m[23], h[23], P[23], h[6], P[6], DT_int_F, DT_F, N_hxrs)
    m[6] = update_mf(m6_req, m[6])
    m[26] = m[25]+m[6]
    # T[24] = temperature( P = P[24], H = h[24])
    f_ip6 = m[6]/m[5]      #IPT1 outlet fraction
    
    #HX E Closed FWH --------------------------------
    print("..HXE", end='')
    # T[22] = temperature( P = P[22], H = h[22])
    # T[26] = temperature( P = P[26], H = h[26])
    m7_req, h[23], h[27], T_H_E[1:3*N_hxrs+1+1], T_C_E[1:3*N_hxrs+1+1] = fwh( m[26], h[26], P[26], m[22], h[22], P[22], h[7], P[7], DT_int_E, DT_E, N_hxrs)
    m[7] = update_mf(m7_req, m[7])
    m[27] = m[26]+m[7]
    # T[23] = temperature( P = P[23], H = h[23])

    m[9] = m[5] - m[6] - m[7] - m[8]    #IPT to LPT
    m[14] = m[9] - m[10] - m[11] - m[12] - m[13]    #LPT to condenser

    f_ip7 = m[7]/(m[5] - m[6])    #IPT2 outlet fraction
    f_ip8 = m[8]/(m[5] - m[6] - m[7])    #IPT3 outlet fraction
    f_lp10 = m[10]/ m[9]    #LPT1 outlet fraction
    f_lp11 = m[11]/(m[9] - m[10])    #LPT2 outlet fraction
    f_lp12 = m[12]/(m[9] - m[10] - m[11])    #LPT3 outlet fraction
    f_lp13 = m[13]/(m[9] - m[10] - m[11] - m[12])    #LPT4 outlet fraction

    # <<<< MJW update estimate of T[1]
    h1old = h[1]
    T1old = T[1]
    h1new = (h[24]*m[24] + h[4]*m[4] - h[25]*m[25])/m[1]
    h[1] = update_mf(h1new, h1old)
    # h[1] = h1old + (h1new-h1old)*0.5

    T[1] = temperature( P = P[1], H = h[1])
    err_iter = abs((T[1] - T1old)/T1old)
    # err_iter = abs((h1new - h1old)/h1old)
    
    # Compute all other temperature states
    for i in range(2,ns):
        T[i] = temperature(H=h[i], P=P[i])

    #Efficiency ===============
    
    W_gen = W_HPT+W_IPT1+W_IPT2+W_LPT1+W_LPT2+W_LPT3+W_LPT4+W_LPT5  #Generator work
    eta_cycle = W_gen/Q_dot_LFR    #Cycle efficiency
    eta_test = (W_gen + Q_dot_cond)/(Q_dot_LFR+W_dot_HP+W_dot_LP)   #Total cycle conservation of energy

    print("    >> Time: {:5.1f}s | Error: {:.{}f} | Power: {:.1f} | eta: {:.4f} | ConsE: {:.4f}".format(time.time()-tstart, err_iter, abs(int(ceil(log10(tol))))+1, W_gen, eta_cycle, eta_test))

    fout.write(log_state(T,P,h,x,m))

    it += 1

fout.close()