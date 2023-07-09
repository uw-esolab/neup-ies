import CoolProp.CoolProp as CP
from state import *
import numpy as np

def convert_C_to_K(temp):
    return temp + 273.15
#######################################################################################################

# -----------------------------------------------------------------
# Example - EES library
# A steam power plant operates on the ideal reheat Rankine cycle.  Steam enters the high pressure turbine
#  at 8 MPa and 500 C and leaves at 3 MPa.  Steam is then reheated at constant pressure to 500 C before 
# it expands to 20 kPa in the low pressure turbine.  Determine the turbine work output, in kJ/kg, and the 
# thermal efficiency of the cycle.  Also show the cycle on a T-s diagram with respect to the saturation lines.
 
# Let's modify this problem to include the effects of the turbine and pump efficiencies and also show the 
# effects of reheat on the steam quality at the low pressure turbine exit.
# -----------------------------------------------------------------

# STATE 0: leaving condenser to pump
# STATE 1: leaving pump to boiler
# STATE 2: leaving boiler to HP turbine
# STATE 3: leaving HP turbine to reheater
# STATE 4: leaving reheater to LP turbine
# STATE 5: leaving LP turbine to condenser


## UnitSystem SI K Pa J

# initialize states as zero
# no. of states
N = 6
# intialize state parameters
P = np.zeros(N)
T = np.zeros(N)
h = np.zeros(N)
s = np.zeros(N)
v = np.zeros(N)
x = np.zeros(N)
states = {'P':P,'T':T,'h':h,'s':s,'v':v,'x':x}

# get pump and turbine isentropic efficiencies
eta_t = 1        # turbine isentropic efficiency
eta_p = 1           # pump isentropic efficiency

# fix known states based on cycle

# known pressures (assuming negligible pressure drop across heat exchangers)
P[2] = 8e6          # known
P[3] = 3e6          # known
P[5] = 20000        # known
P[1] = P[2]
P[0] = P[5]
P[4] = P[3]

# known temperatures
T[2] = convert_C_to_K(500)
T[4] = convert_C_to_K(500)

# known quality
x[0] = 0            # quality of steam out of condenser

#--------------------------------------------------------------------------------
# Pump analysis (0-->1)

# get state 0
h[0] = enthalpy(X=x[0],P=P[0])
v[0] = volume(H=h[0],P=P[0])
s[0] = entropy(H=h[0],P=P[0])
T[0] = temperature(H=h[0],P=P[0])
# get pump power required
W_pump_s = v[0] * (P[1] - P[0])     # isentropic pump work assuming constant specific volume
W_pump = W_pump_s/eta_p
# get state 1
h[1] = h[0] + W_pump
v[1] = volume(H=h[1],P=P[1])
s[1] = entropy(H=h[1],P=P[1])
T[1] = temperature(H=h[1],P=P[1])

#---------------------------------------------------------------------------------
# Boiler Analysis (1-->2)

# get state 2
h[2] = enthalpy(T=T[2],P=P[2])
v[2] = volume(H=h[2],P=P[2])
s[2] = entropy(T=T[2],P=P[2])

# get boiler heat in
Q_in_boiler = h[2] - h[1]

#---------------------------------------------------------------------------------
# HP turbine analysis (2-->3)

# get HP turbine power produced
ss_3 = s[2]
hs_3 = enthalpy(S=ss_3,P=P[3])

W_HPT_s = h[2] - hs_3
W_HPT = W_HPT_s * eta_t

# get state 3
h[3] = h[2] - W_HPT
T[3] = temperature(P=P[3],H=h[3])
s[3] = entropy(P=P[3],H=h[3])
v[3] = volume(P=P[3],H=h[3])

#---------------------------------------------------------------------------------
# Reheater (3-->4)

# get state 4
h[4] = enthalpy(T=T[4],P=P[4])
v[4] = volume(H=h[4],P=P[4])
s[4] = entropy(H=h[4],P=P[4])
T[4] = temperature(H=h[4],P=P[4])

# get Q_in_reheat
Q_in_reheat = h[4] - h[3]

#---------------------------------------------------------------------------------
# LP turbine analysis (4-->5)

# get LPT power produced
ss_5 = s[4]
hs_5 = enthalpy(S=ss_5,P=P[5])

W_LPT_s = h[4] - hs_5
W_LPT = W_LPT_s * eta_t

# get state 5
h[5] = h[4] - W_LPT
x[5]=quality(H=h[5],P=P[5])
T[5]=temperature(H=h[5],P=P[5])
s[5]=entropy(H=h[5],P=P[5])
#---------------------------------------------------------------------------------
# Condenser (5-->0)
Q_out = h[5] - h[0]

#---------------------------------------------------------------------------------
# Cycle Statistics
W_net = W_HPT + W_LPT - W_pump
eff = W_net/(Q_in_boiler + Q_in_reheat)
