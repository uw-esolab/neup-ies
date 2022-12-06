import CoolProp.CoolProp as CP
from state import *
from components import *
import numpy as np
import time
from math import log10,ceil
#######################################################################################################

# -----------------------------------------------------------------
# Example - Implicit Rankine cycel with Feedwater heater
# A steam power plant operates on the ideal reheat Rankine cycle. Steam enters the high pressure turbine
# at 8 MPa and 500 C and leaves at 3 MPa. The drain from the turbine is sent to a high pressure and low 
# pressure feed water heater cascaded backwards to the condenser. Determine the turbine work output, in kJ/kg, and the 
# thermal efficiency of the cycle.  Also show the cycle on a T-s diagram with respect to the saturation lines.
 
# Let's modify this problem to include the effects of the turbine and pump efficiencies
# -----------------------------------------------------------------

# STATE 1: Boiler to Turbine
# STATE 2: Turbine to FWH_A
# STATE 3: Turbine to FWH_B
# STATE 4: Turbine exit to Condenser
# STATE 5: Condenser to Pump
# STATE 6: Pump to FWH_B (main flow)
# STATE 7: FWH_B ro FWH_A (main flow)
# STATE 8: FWH_A to boiler
# STATE 9: FWH_B to condenser (drain flow)
# STATE 10: FWH_A to FWH_B (drain flow)

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

# ---------------------------------------------------------------------------------
## INPUTS FOR SOLVER
DT = 5.5 #[K]      #Target Approach Temperature 
DT_int = 3 #[K]
N_hxrs = 5 #[-]    #Number of heat exchangers
ns = 10  #number of states used + 1
T_r = 25+273.15    #condenser Reservoir Temperature
# ---------------------------------------------------------------------------------

# start writing log to file 
fout = open('iter-log.csv','w')
fout.write(log_header(ns, "T P h x m".split(' ')))

# initialize arrays for cycle
m = np.zeros(ns)                        # mass flow
T_H_A = np.zeros(3*N_hxrs + 2)          # Temperatures array for FWH_A hot side
T_C_A = np.zeros(3*N_hxrs + 2)          # Temperatures array for FWH_A cold side
T_H_B = np.zeros(3*N_hxrs + 2)          # Temperatures array for FWH_B hot side
T_C_B = np.zeros(3*N_hxrs + 2)          # Temperatures array for FWH_A cold side
h = np.zeros(ns)                        # enthalpy
P = np.zeros(ns)                        # pressure
T = np.zeros(ns)                        # temperature
m_t = np.zeros(ns)                      # intermediate mass flow during turbien drain
x = np.zeros(ns)                        # quality

loc = np.zeros(3*N_hxrs+2)              # location for FWH temp distribution
for i in range(1,(3*N_hxrs+1)+1):
    loc[i] = (i-1)/(3*N_hxrs)

PR = 2.657 #[-] even pressure across turbines assumption
# Mass Flow Rates ==================================================================

m[1] = 477.2 #[kg/s]

# initial calcualtions based on guess values
f_A = 0.001                             # fraction to FWH_A
f_B = 0.1                               # fraction to FWH_B
m[6] = m[1]                             # intial guess for mass flow to FWH_B
h[10] = 436116
h[6] = 663252

P[1] = 30e6                         # Boiler outlet pressure
T[1] = 632 + 273.15 #[K]            # Boiler outlet temperature
h[1] = enthalpy(T=T[1],P=P[1])      # Boiler outlet enthalpy

eta_T1 = 0.90 #[-]                  # T1 isentropic efficiency
eta_T2 = 0.90 #[-]                  # T2 isentropic efficiency
eta_T3 = 0.90 #[-]                  # T2 isentropic efficiency
eta_P = 0.85 #[-]                   # Pump efficiency

DT_A = 5.5 #[K]                     # FWH_A approach temperature difference
DT_int_A = 3 #[K]                   # FWH_A pinch point temperature difference
DT_B = 5.5 #[K]                     # FWH_B approach temperature difference
DT_int_B = 3 #[K]                   # FWH_B pinch point temperature difference
T_r = 25 + 273.15 #[K]              # Condenser Reservoir Temperature

P[6] = 17.24 * 1e5                  # Pump outlet pressure
x[4] = 0.9 #[-]                     # Turbine outlet quality
P[4] = 5000 #[Pa]                   # Turbine outlet pressure
P[3] = P[4] * PR                    # Turbine FWH_B outlet pressure
P[2] = P[3] * PR                    # Turbien FWH_A outlet pressure
# Initial guess for T[8]
T[8] = 611.4 #[C]

it = 0
it_max = 50
err_iter = 999.
tol = 0.0001

while((err_iter > tol and it<it_max) or it<3 ):

    print("Iter {:d}: ".format(it), end='')

    # mass flows known -----------------------------------------------------------------
    m[8] = m[1]                         # Boiler
    m[2] = f_A * m[1]                   # FWH_A from turbine
    m[3] = (m[1] - m[2]) * f_B          # FWH_B from turbine
    m[4] = m[1] - m[2] - m[3]           # From turbine to condenser
    #-----------------------------------------------------------------------------------

    #LFR Boiler -----------------------------------------
    delta_P_boiler = 5*286107  #[Pa]    #Assumption of pressure drop accross boiler (Jacob Wenner)
    P[8] = P[1] + delta_P_boiler    #Boiler inlet pressure
    
    # {T[1] = 340+273.15    #Boiler inlet temperature}
    h[8] = enthalpy( T = T[8], P = P[8])    #Boiler inlet enthalpy
    
    Q_dot_LFR = 950e6    #Heat transfer from LFR
    
    h[1] = (m[8]*h[8] + Q_dot_LFR)/m[1]    #Energy balance on Boiler

    # Constants for components (isentropic efficiency and approach temperature differences)
    # Once again ask if move outside loop

    # Turbine ---------------------------------------------

    # T1
    h[2], W_T1 = turbine( m[1], h[1], P[1], P[2], eta_T1)
    T[2] = temperature(H=h[2],P=P[2])

    # T2
    m_t[2] = m[1] - m[2]                # turbine flow from 2 to 3
    h[3], W_T2 = turbine( m_t[2], h[2], P[2], P[3], eta_T2)
    T[3] = temperature(H=h[3],P=P[3])

    # T3
    m_t[3] = m_t[2] - m[3]              # turbine flow from 3 to 4
    h[4], W_T3 = turbine( m_t[3], h[3], P[3], P[4], eta_T3)
    T[4] = temperature(H=h[4],P=P[4])

    # get mass flows and pressure from equation
    m[10] = m[2]
    P[10] = P[2]
    P[7] = P[6]

    
    # HX B closed FWH -------------------------------------
    print("..HXB", end='')
    m3_req, h[7], h[9], T_H_B[1:3*N_hxrs+1+1], T_C_B[1:3*N_hxrs+1+1] = fwh(m[10], h[10], P[10], m[6], h[6], P[6], h[3], P[3], DT_int, DT, N_hxrs)
    m[3] = update_mf(m7_req, m[7])

    # HX A closed FWH -------------------------------------
    h8_old = h[8]
    m2_req, h[8], h[10], T_H_A[1:3*N_hxrs+1+1], T_C_A[1:3*N_hxrs+1+1] = fwh(0, 0, 0, m[7], h[7], P[7], h[2], P[2], DT_int, DT, N_hxrs)
    m[2] = update_mf(m2_req, m[2])

    # Condenser -------------------------------------------

    # add energy balance into condenser
    



    # Pump ------------------------------------------------

    # Add pump energy balance



    # update fractions ------------------------------------






    # <<<< MJW update estimate of T[1]
    h8old = h[8]
    T8old = T[8]
    h8new = (h[7]*m[7] + h[2]*m[2] - h[10]*m[10])/m[8]
    h[8] = update_mf(h8new, h8old)
    
    T[8] = temperature( P = P[8], H = h[8])
    err_iter = abs((T[8] - T8old)/T8old)
    # err_iter = abs((h8new - h8old)/h8old)
    
    # Compute all other temperature states
    for i in range(1,ns):
        T[i] = temperature(H=h[i], P=P[i])

    # Effieciency ==========================================

    #W_gen = W_T1 + W_T2                 # Generator work
    #eta_cycle = W_gen / Q_dot_LFR       # Cycle efficiency
    #eta_test = (W_gen + Q_dot_cond) / (Q_dot_LFR = W_dot_pump)

    #print("    >> Time: {:5.1f}s | Error: {:.{}f} | Power: {:.1f} | eta: {:.4f} | ConsE: {:.4f}".format(time.time()-tstart, err_iter, abs(int(ceil(log10(tol))))+1, W_gen, eta_cycle, eta_test))

    fout.write(log_state(T,P,h,x,m))

    it += 1

fout.close()