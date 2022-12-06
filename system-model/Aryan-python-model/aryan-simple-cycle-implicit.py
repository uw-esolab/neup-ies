from distutils.log import error
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

# out of turbine --> 90% quality 5 mbar
# extremely simple model
# -----------------------------------------------------------------

# states
# 1: going to boiler
# 2: boiler to turbine
# 3: turbine to feedwater drain
# 4: turbine to condenser
# 5: mixed to condenser
# 6: condenser to pump
# 7: pump to fwh (main flow)
# 8: fwh(drain) to mixer

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

# solver inputs
DT = 5.5 #[K]
ns = 9             # number of states + 1
T_r = 25+273.15    #condenser Reservoir Temperature

# start writing log to file 
fout = open('iter-log.csv','w')
fout.write(log_header(ns, "T P h x m".split(' ')))

# initialize arrays for cycle
m = np.zeros(ns)                        # mass flow
h = np.zeros(ns)                        # enthalpy
P = np.zeros(ns)                        # pressure
T = np.zeros(ns)                        # temperature
m_t = np.zeros(ns)                      # intermediate mass flow during turbien drain
x = np.zeros(ns)                        # quality
s = np.zeros(ns)                        # entropy

m[1] = 477.2 #[kg/s]                    # solve using state points
Q_dot_LFR = 950e6                       # known value

# intial guess values
f = 0.1

# get P_sat 
# replace h[8] with the DT value. DT value gives T[8]. from T[8] and P[8] = P[3] --> h[8]
h[8] = 863252

# known state values
T[1] = 611.4 #[K]
P[2] = 30e6 #[Pa]                       # Boiler outlet pressure
T[2] = 632 + 273.15                     # Boiler outlet temperature
h[2] = enthalpy(T=T[2],P=P[2]) 
delta_P_boiler = 5*286107 #[Pa]         # Assumption of pressure drop across the boiler

eta_T = 0.90 #[-]
eta_P = 0.85 #[-]

x[4] = 0.9 #[-]                         # Turbine outlet quality
P[4] = 5000 #[Pa]                       # Turbine outlet presssure
P[7] = 17.24 * 1e5                  # Pump outlet pressure
PR = sqrt(P[2]/P[4]) #[-]                         # updated turbine pressure ratio
P[3] = PR * P[4]

# intial guess value
it = 0
it_max = 50
err_iter = 999.
tol = 0.0001

## start cycle solving
while((err_iter > tol and it<it_max) or it<3 ):

    # ----------------------------------------------------------------------------------------------    
    # boiler

    P[1] = P[2] + delta_P_boiler #[Pa]
    h[1] = enthalpy(P=P[1],T=T[1])
    #P[7] = P[1]
    #Q_dot_LFR = m[2]*h[2] - m[1]*h[1]
    # update to solve for mass flow
    m[1] = Q_dot_LFR/(h[2]-h[1])
    m[2] = m[1]
    # get mass flows
    m[3] = f * m[2]

    # ----------------------------------------------

    # turbine
    # T1
    h[3], W_T1 = turbine( m[2], h[2], P[2], P[3], eta_T)
    T[3] = temperature(H=h[3],P=P[3])

    # T2
    m_t[2] = m[2] - m[3]
    m[4] = (1-f)*m[2]
    h[4], W_T2 = turbine( m_t[2], h[3], P[3], P[4], eta_T)
    T[4] = temperature(H=h[4],P=P[4])
    s[4] = entropy(H=h[4],P=P[4])
    # calculate quality at 4

    # get mass flows and pressure

    # get h[8]
    m[8] = m[3]
    P[8] = P[3]
    s[8] = entropy(P=P[8],H=h[8])

    # mixing
    # get h[8] from approach temperature difference
    
    #----#

    # enthalpy balance
    m[5] = m[4] + m[8]
    h[5] = (m[4]*h[4] + m[8]*h[8])/(m[5])
    # assuming isobaric mixing
    P[5] = P[4]

    # condenser
    x[6] = 0
    P[6] = P[5]
    # switch to energy balance. Get h[6] from quality and pressure
    #h[6] = enthalpy(P=P[6],x=x[6])
    h[6], P[6], Q_dot_cond, epsilon_cond = condenser( m[5], h[5], P[5], T_r, DT)
    m[6] = m[5]


    # pump
    h[7], W_dot_HP = pump( m[6], h[6], P[6], P[7], eta_P)

    # fwh

    # do full energy on fwh
    # 
    h1old = h[1]
    T1old = T[1]
    h1new = (h[3]*m[3] + h[7]*m[7] - h[8]*m[8])/m[1]

    # update fraction of 
    h[1] = update_mf(h1new, h1old)
    # h[1] = h1old + (h1new-h1old)*0.5

    T[1] = temperature( P = P[1], H = h[1])
    err_iter = abs((T[1] - T1old)/T1old)

    # update fractions
    f = f + 0.1*err_iter
    for i in range(2,ns):
            T[i] = temperature(H=h[i], P=P[i])

    #Efficiency ===============
        
    W_gen = (W_T1 + W_T2 - W_dot_HP)
    eta_cycle = W_gen/Q_dot_LFR    #Cycle efficiency
    eta_test = (W_gen + Q_dot_cond)/(Q_dot_LFR+W_dot_HP)   #Total cycle conservation of energy

    print("    >> Time: {:5.1f}s | Error: {:.{}f} | Power: {:.1f} | eta: {:.4f} | ConsE: {:.4f}".format(time.time()-tstart, err_iter, abs(int(ceil(log10(tol))))+1, W_gen, eta_cycle, eta_test))

    fout.write(log_state(T,P,h,x,m))

    it += 1

#fout.close()

end_state = {'T':[T[7],T[1],T[3],T[8]],'P':[P[7],P[1],P[3],P[8]],'m':[m[1],m[3]]}
print(end_state)
print(W_dot_HP)
print(W_T2)