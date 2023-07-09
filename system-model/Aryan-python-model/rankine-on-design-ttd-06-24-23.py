from state import *
from components_v2 import *
import numpy as N
import time
from math import *

##  OFF Design System Model for Supercritical Rankine Steam Cycle =================================
##  Based off of Brian's EES model --> rankine-on-design ==========================================

# Define function for logging and reporting data
# --------------------------
def log_header(N,names):
    """Writes the headers for the state output file."""
    hh = []
    for name in names:
        h = []
        for i in range(1,N):
            h.append("{:s}{:d}".format(name,i))
        hh.append(h)
    return ','.join([','.join(h) for h in hh]) + '\n'

def log_state(*args):
    """Returns the state properties for the converged system."""
    # T,P,s,h,x,m
    return ','.join([','.join([str(v) for v in arr[1:]]) for arr in args]) + '\n'

step = lambda xold,xnew: xold + (xnew-xold)*0.66
# --------------------------

# Define function for updating mass flow with every iteration
def update_mf(m_dot_old, err):
    # TODO: update convergence parameters based on CoE.
    gr = (1 + sqrt(5))/2
    return m_dot_old * ()
tstart = time.time()

# Model start

# Select mode of operation --> O = OFF, C = 100% Charging, D = 100% Discharging
# mode = input("Please select operation mode:")
mode = "O"                                                         

#-----------------------------------------------#
# TODO: update values
DT = 5.55 #[K]                                                      # Target Approach Temperature 
N_hxrs = 20 #[-]                                                    # Number of sub heat exchangers
ns = 33                                                             # number of states used + 1
#-----------------------------------------------#
# Start writing to iteration log file 
fout = open('iter-log2.csv','w')
fout.write(log_header(ns, "T P h x m".split(' ')))                  # provide headers to log file

# Initialize arrays
m = N.zeros(ns)                                                     # mass flow array
m_t = N.zeros(ns)                                                   # turbine mass flow array
m_s = N.zeros(ns)                                                   # salt mass flow array
T_H_A = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX A
T_C_A = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX A
T_H_B = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX B
T_C_B = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX B
T_H_C = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX C
T_C_C = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX C
T_H_D = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX D
T_C_D = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX D
T_H_E = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX E
T_C_E = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX E
T_H_F = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX F
T_C_F = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX F
T_H_G = N.zeros(3*N_hxrs+2)                                         # hot side temperature at HX G
T_C_G = N.zeros(3*N_hxrs+2)                                         # cold side temperature at HX G
h = N.zeros(ns)  * float('NaN')                                                     # enthalpy array
P = N.zeros(ns)                                                     # pressure array
T = N.zeros(ns)                                                     # temperature array
x = N.zeros(ns)                                                     # quality array

loc = N.zeros(N_hxrs*3+2)                                           # position array across each hx
# Populate location array based on number of sub heat exchangers
for i in range(1,(3*N_hxrs+1)+1):
    loc[i] = (i-1)/(3*N_hxrs)

# Initialize heat exchanger parameters
if (mode.lower() == "o" or mode.lower() == "off"):
    # all temperatures in Kelvin
    TTD_A = 2.72                                                    # terminal temperature difference
    DT_HX_A = 5.07                                                  # approach temperature difference
    TTD_B = 2.72
    DT_HX_B = 5.03
    TTD_C = 2.74
    DT_HX_C = 5.01
    TTD_D = 2.90
    DT_HX_D = 4.99
    TTD_E = 2.43
    DT_HX_E = 5.05
    TTD_F = 2.28
    DT_HX_F = 5.00
    DT_HX_G = 4.37
elif (mode.lower() == "c" or mode.lower() == "charging"):
    # all temperatures in Kelvin
    TTD_A = 2.63                                                    # terminal temperature difference
    DT_HX_A = 4.93                                                  # approach temperature difference
    TTD_B = 2.62
    DT_HX_B = 4.89
    TTD_C = 2.64
    DT_HX_C = 4.86
    TTD_D = 2.79
    DT_HX_D = 4.86
    TTD_E = 2.86
    DT_HX_E = 6.57
    TTD_F = 2.84
    DT_HX_F = 6.83
    DT_HX_G = 6.40
    m_s[7] = 190 #[kg/s]                                            # salt mass flow through W2S
elif (mode.lower == "d" or mode.lower() == "discharging"):
    # all temperatures in Kelvin
    TTD_A = 3.00                                                    # terminal temperature difference
    DT_HX_A = 5.55                                                  # approach temperature difference
    TTD_B = 3.00
    DT_HX_B = 5.55
    TTD_C = 3
    DT_HX_C = 5.55 
    TTD_D = 3
    DT_HX_D = 5.55 
    TTD_E = 3
    DT_HX_E = 5.55 
    TTD_F = 3
    DT_HX_F = 5.55 
    DT_HX_G = 5.55
    m_s[4] = 525 #[kg/s]                                            # salt mass flow through S2W
    m[33] = 100 #[kg/s]                                             # steam mass flow through W2S
else:
    raise Exception("Operation Mode not recognized. Must be either 'O', 'C', or 'D'.")

# Component Efficiencies
eta_HPT = 0.90  #[-]                                                # HPT isentropic efficiency
eta_IPT1 = 0.90 #[-]                                                # IPT1 isentropic efficiency
eta_IPT2 = 0.90 #[-]                                                # IPT2 isentropic efficiency
eta_IPT3 = 0.90 #[-]                                                # IPT3 isentropic efficiency
eta_LPT1 = 0.90 #[-]                                                # LPT1 isentropic efficiency
eta_LPT2 = 0.90 #[-]                                                # LPT2 isentropic efficiency
eta_LPT3 = 0.90 #[-]                                                # LPT3 isentropic efficiency
eta_LPT4 = 0.90 #[-]                                                # LPT4 isentropic efficiency
eta_LPT5 = 0.90 #[-]                                                # LPT5 isentropic efficiency
eta_LP = 0.90   #[-]                                                # LP pump efficiency

# initialize pressure out of LP turbine

"""TODO: Start by defining all knowns. Then add in component_v2 library. 
Figure out what how many guess values needed. --> Should be 1. Otherwise need more implicit functions. 
Refer to Brian's code and setup model."""

# Knowns
P[14] = 5000 #[Pa]                                                 # Inlet pressure of LP Turbine
P[16] = 2.544e6 #[Pa]                                               # Outlet pressure of LP Turbine

# Get Mass Flow Rates ==========
 
# LFR Boiler
P[2] = 3.3e7                                                        # Boiler outlet pressure
T[2] = 632+273.15                                                   # Boiler outlet temperature
h[2] = enthalpy( T = T[2], P = P[2])                                # Boiler outlet enthalpy
P[1] = P[2]                                                         # Boiler inlet pressure, assuming no pressure loss across boiler
T[1] = 340+273.15                                                   # Boiler inlet temperature --> KNOWN
h[1] = enthalpy( T = T[1], P = P[1])                                # Boiler inlet enthalpy
Q_dot_LFR = 950e6                                                   # Heat transfer from LFR

m[1] = Q_dot_LFR/(h[2]-h[1])                                        # Get cycle mass flow rate
m[2] = m[1]
h[32] = h[31]
m[21] = m[1]

# High Pressure Turbines

# HPT 1
P[3] = 2.365e7 #[Pa]                                                # 

# Initial guess for m[15]
m[15] = 300 #[kg/s]

it = 0
it_max = 50
err_iter = 999.
tol = 0.0001

while((err_iter > tol and it<it_max) or it<3 ):

    print("Iter {:d}: ".format(it), end='')

    # Temperature Estimation out of LP Pump
    P[15] = P[14]                                                       # No Pressure Drop across Condenser
    T[15] = temperature(P=P[15],X=0)                                    # LP Pump inlet temperature
    h[15] = enthalpy(X=0, P=P[15])                                   # LP Pump inlet enthalpy
    h[16], W_dot_LP = pump( m[15], h[15], P[15], P[16], eta_LP)
    T[16] = temperature(P=P[16],H=h[16])                                # LP Pump outlet temperature
    m[16] = m[15]
    # Get FWH temperature rises assuming even temp rise across FWH
    n_fwh = 8                                                           # Number of feed water heaters
    T_rise = (T[1] - T[16])/n_fwh
    T[17] = T[16] + T_rise
    T[18] = T[17] + T_rise
    T[19] = T[18] + T_rise
    T[20] = T[19] + T_rise
    T[22] = T[20] + T_rise
    T[23] = T[22] + T_rise
    T[24] = T[23] + T_rise
    T[32] = T[1]
    # Pressure in FW
    P[17] = P[16]                                                       # Low Pressure FW
    P[18] = P[17]
    P[19] = P[18]
    P[20] = P[19]

    P[22] = P[1]                                                        # High Pressure FW
    P[23] = P[22]
    P[24] = P[23]
    P[32] = P[24]

    # Heat Transfer required in hp FW
    for i in range(18,25):
        if i != 21:
            h[i] = enthalpy(T=T[i],P=P[i])
    h[32] = enthalpy(P=P[32],T=T[32])

    # HP Turbines
    m[17] = m[16]
    m[18] = m[17]
    m[19] = m[18]
    m[20] = m[19]
    m[22] = m[21]
    m[23] = m[22]
    m[24] = m[23]

    Q_dot_HXE = m[22]*(h[23]-h[22])
    Q_dot_HXF = m[23]*(h[24]-h[23])
    Q_dot_HXG = m[24]*(h[32]-h[24])

    P[6] = pressure(T = (T[24]+TTD_F), X = 0.1)
    P[5] = (P[6]+P[3])/2
    P[7] = pressure(T = (T[23]+TTD_E), X = 0.1)
    P[8] = pressure(T = T[22], X = 0.1)
    P[9] = P[8]
    P[10] = pressure(T = (T[20]+TTD_D), X = 0.1)
    P[11] = pressure(T = (T[19]+TTD_C), X = 0.1)
    P[12] = pressure(T = (T[18]+TTD_B), X = 0.1)
    P[13] = pressure(T = (T[17]+TTD_A), X = 0.1)

    # HP energy balance for turbine extraction mass flows

    # HPT 1
    # Define HPT1 states for OFF design
    P[4] = P[3]
    P[25] = P[4]
    T[25] = T[24]+DT_HX_G
    h[25] = enthalpy(T = T[25], P = P[25])
    h[4], h_HPT1, m[4], m_HPT1, x_HPT1, W_dot_HPT1 = turbine(m[2], h[2], P[2], P[3], eta_HPT, Q_dot_HXG, h[25], 0, 0)
    Q_dot_W2S = 0
    Q_dot_S2W = 0
    m[3] = m[4]
    T[4] = T[3]
    m[25] = m[4]
   
    h[3] = h[4]
    T[3] = temperature(P = P[3], H = h[3])
    x[3] = quality(P = P[3], H = h[3])
    print(x[3])

    # HPT 2
    h[5], h_HPT2, m[5], m_HPT2, x_HPT2, W_dot_HPT2 = turbine(m_HPT1, h_HPT1, P[3], P[5], eta_HPT, 0,0,0,0)
    T[5] = temperature(P = P[5], H = h[5])
    x[5] = quality(P = P[5], H = h[5])
    m[5] = m_HPT2
    h[5] = h_HPT2
    
    # IPT 1
    T[26] = T[23]+DT_HX_F
    P[26] = P[6]
    h[26] = enthalpy(T = T[26], P = P[26])
    h[6], h_IPT1, m[6], m_IPT1, x_IPT1, W_dot_IPT1 = turbine(m[5], h[5], P[5], P[6], eta_IPT1, Q_dot_HXF, h[26], m[25], h[25])
    m[26] = m[25]+m[6]
    T[6] = temperature(H = h[6], P = P[6])
    x[6] = quality(P = P[6], H = h[6])
 
    # IPT 2
    T[27] = T[22]+DT_HX_E
    P[27] = P[7]
    h[27] = enthalpy(T = T[27], P = P[27])
    h[7], h_IPT2, m[7], m_IPT2, x_IPT2, W_dot_IPT2 = turbine(m_IPT1, h_IPT1, P[6], P[7], eta_IPT2, Q_dot_HXE, h[27], m[26], h[26])
    m[27] = m[26]+m[7]
    T[7] = temperature(P = P[7], H = h[7])
    x[7] = quality(P = P[7], H = h[7])
 
    # IPT3
    P[21] = min(P[8],P[20],P[27])
    m[21] = m[1]
    m[8] = -(m[20] + m[27]) + m[21]
    h[8], W_dot_IPT3 = turbine_simple(m_IPT2, h[7], P[7], P[8], eta_HPT)
    h[21] = (m[20]*h[20] + (m[8]*h[8] + m[27]*h[27]))/m[21]
    h[22], W_dot_HP = pump(m[21], h[21], P[21], P[22], eta_IPT3)
    T[21] = temperature(P = P[21], H = h[21])
    T[8] = temperature(P = P[8], H = h[8])
    x[8] = quality(P = P[8], H = h[8])
    x_IPT3 = x[8]
 
    # LP Turbines     
    h[17] = enthalpy(P=P[17],T=T[17])
    Q_dot_HXA = m[16]*(h[17]-h[16])
    Q_dot_HXB = m[17]*(h[18]-h[17])
    Q_dot_HXC = m[18]*(h[19]-h[18])
    Q_dot_HXD = m[19]*(h[20]-h[19])

    # LPT 1
    h[9] = h[8]
    T[9] = T[8]
    m[9] = m[5] - m[6] - m[7] - m[8]
    T[28] = T[19] + DT_HX_D
    P[28] = P[10]
    h[28] = enthalpy(T = T[28], P = P[28])
    h[10], h_LPT1, m[10], m_LPT1, x_LPT1, W_dot_LPT1 = turbine(m[9], h[9], P[9], P[10], eta_LPT1, Q_dot_HXD, h[28], 0, 0)
    T[10] = temperature(P = P[10], H = h[10])
    x[10] = quality(P = P[10], H = h[10])
    m[28] = m[10]
 
    # LPT 2
    T[29] = T[18]+DT_HX_C
    P[29] = P[11]
    h[29] = enthalpy(T = T[29], P = P[29])
    print(m_LPT1)
    h[11], h_LPT2, m[11], m_LPT2, x_LPT2, W_dot_LPT2 = turbine(m_LPT1, h_LPT1, P[10], P[11], eta_LPT2, Q_dot_HXC, h[29], m[28], h[28])
    T[11] = temperature(P = P[11], H = h[11])
    x[11] = quality(P = P[11], H = h[11])
    m[29] = m[28]+m[11]
 
    #LPT 3
    T[30] = T[17]+DT_HX_B
    P[30] = P[12]
    h[30] = enthalpy(T = T[30], P = P[30])
    h[12], h_LPT3, m[12], m_LPT3, x_LPT3, W_dot_LPT3 = turbine(m_LPT2, h_LPT2, P[11], P[12], eta_LPT3, Q_dot_HXB, h[30], m[29], h[29])
    T[12] = temperature(P = P[12], H = h[12])
    x[12] = quality(P = P[12], H = h[12])
    m[30] = m[29]+m[12]
 
    # LPT 4
    T[31] = T[16]+DT_HX_A
    P[31] = P[13]
    h[31] = enthalpy(T = T[31], P = P[31])
    h[13], h_LPT4, m[13], m_LPT4, x_LPT4, W_dot_LPT4 = turbine(m_LPT3, h_LPT3, P[12], P[13], eta_LPT4, Q_dot_HXA, h[31], m[30], h[30])
    T[13] = temperature(P = P[13], H = h[13])
    x[13] = quality(P = P[13], H = h[13])
    m[31] = m[30]+m[13]
  
    # LPT 5
    h[14],W_dot_LPT5 = turbine_simple(m_LPT4, h_LPT4, P[13], P[14], eta_LPT5)
    T[14] = temperature(P = P[14], H = h[14])
    x[14] = quality(P = P[14], H = h[14])
    m[14] = m[9]-m[10]-m[11]-m[12]-m[13]
    m[15] = m[14]+m[31]
 
    # condenser
    h_cond_in = (h[14]*m[14] + h[31]*m[31]) / (m[14]+m[31])
    Q_dot_cond = (m[14]+m[31])*(h_cond_in-h[15])

    # Efficiencies
    W_dot_HPT = W_dot_HPT1+W_dot_HPT2
    W_dot_IPT = W_dot_IPT1+W_dot_IPT2+W_dot_IPT3
    W_dot_LPT = W_dot_LPT1+W_dot_LPT2+W_dot_LPT3+W_dot_LPT4+W_dot_LPT5
    W_dot_gen = W_dot_HPT1+W_dot_HPT2+W_dot_IPT1+W_dot_IPT2+W_dot_IPT3+W_dot_LPT1+W_dot_LPT2+W_dot_LPT3+W_dot_LPT4+W_dot_LPT5
    # Conservation of energy check
    eta_test = (W_dot_gen+Q_dot_cond+Q_dot_W2S)/(Q_dot_LFR+W_dot_LP+W_dot_HP+Q_dot_S2W)

    # <<<< MJW update estimate of T[1]
    m15old = m[15]
    m15new = m[14] + m[31]
    m[15] = update_mf(m15new, m15old)
    print(m15new, m15old)
    # h[1] = h1old + (h1new-h1old)*0.5
    err_iter = abs((1 - eta_test)/eta_test)
    # err_iter = abs((h1new - h1old)/h1old)
    
    # Compute all other temperature states
    for i in range(2,ns):
        if (h[i] != 0) and (P[i] != 0):
            T[i] = temperature(H=h[i], P=P[i])

    #Efficiency ===============
    
    eta_cycle = W_dot_gen/Q_dot_LFR    # Cycle efficiency

    print("    >> Time: {:5.1f}s | Error: {:.{}f} | Power: {:.1f} | eta: {:.4f} | ConsE: {:.4f}".format(time.time()-tstart, err_iter, abs(int(ceil(log10(tol))))+1, W_dot_gen, eta_cycle, eta_test))

    fout.write(log_state(T,P,h,x,m))

    it += 1

fout.close()