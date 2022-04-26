from state import *

# ######################################################################################################

# -----------------------------------------------------------------
# Example - EES library
# A steam power plant operates on the ideal reheat Rankine cycle.  Steam enters the high pressure turbine
#  at 8 MPa and 500 C and leaves at 3 MPa.  Steam is then reheated at constant pressure to 500 C before 
# it expands to 20 kPa in the low pressure turbine.  Determine the turbine work output, in kJ/kg, and the 
# thermal efficiency of the cycle.  Also show the cycle on a T-s diagram with respect to the saturation lines.
 
# Let's modify this problem to include the effects of the turbine and pump efficiencies and also show the 
# effects of reheat on the steam quality at the low pressure turbine exit.
# -----------------------------------------------------------------

N = 6
P = [0. for i in range(N)]
T = [0. for i in range(N)]
h = [0. for i in range(N)]
s = [0. for i in range(N)]
v = [0. for i in range(N)]
x = [0. for i in range(N)]

# $UnitSystem SI C kPa kJ

P[5] = 20000    #[Pa]
P[2] = 8e6      #[Pa]
T[2] = 773.15   #[K]
P[3] = 3e6      #[Pa]
T[4] = 773.15   #[K]
Eta_t = 0.86 	#"Turbine isentropic efficiency"
Eta_p = 0.5	#"Pump isentropic efficiency"

 
# "!Pump analysis"
P[0] = P[5]
P[1] = P[2]
x[0]=0 	#"Sat'd liquid"
# This is the first way of getting the state, which only calls the state functions once. It's faster.
state_0 = State("pump", P=P[0], X=x[0])
h[0] = state_0.H.m
v[0] = state_0.v.m 
s[0] = state_0.S.m 
T[0] = state_0.T.m

W_p_s=v[0]*(P[1]-P[0])	#"SSSF isentropic pump work assuming constant specific volume"
W_p=W_p_s/Eta_p  
h[1]=h[0]+W_p   	#"SSSF First Law for the pump"
# This is the second way of getting state properties, which creates a state object each time a function is called. It's slower,
# but more like EES
v[1]=volume(P=P[1],H=h[1])
s[1]=entropy(P=P[1],H=h[1])
T[1]=temperature(P=P[1],H=h[1])
 
# "!High Pressure Turbine analysis"
h[2]=enthalpy(T=T[2],P=P[2])
s[2]=entropy(T=T[2],P=P[2])
v[2]=volume(T=T[2],P=P[2])
s_s_4=s[2]
hs_4=enthalpy(S=s_s_4,P=P[3])
Ts_4=temperature(S=s_s_4,P=P[3])
#"Definition of turbine efficiency"
h[3] = h[2] - (h[2]-hs_4)*Eta_t

T[3]=temperature(P=P[3],H=h[3])
s[3]=entropy(T=T[3],P=P[3])
v[3]=volume(S=s[3],P=P[3])
# "SSSF First Law for the high pressure turbine"
W_t_hp = h[2] - h[3]
 
# "!Low Pressure Turbine analysis"
P[4]=P[3]
s[4]=entropy(T=T[4],P=P[4])
h[4]=enthalpy(T=T[4],P=P[4])
s_s_6=s[4]
hs_6=enthalpy(S=s_s_6,P=P[5])
Ts_6=temperature(S=s_s_6,P=P[5])
vs_6=volume(S=s_s_6,P=P[5])
# "Definition of turbine efficiency"
h[5] = h[4]-Eta_t*(h[4]-hs_6)
# "SSSF First Law for the low pressure turbine"
W_t_lp = h[4] - h[5]	
x[5]=quality(H=h[5],P=P[5])
 
# "!Boiler analysis"
Q_in = h[2]+h[4] - (h[1]+h[3])	#"SSSF First Law for the Boiler"
 
# "!Condenser analysis"
Q_out = h[5] - h[0]	#"SSSF First Law for the Condenser"
T[5]=temperature(H=h[5],P=P[5])
s[5]=entropy(H=h[5],P=P[5])
# x6s$=' ('||phase$(H=h[5],P=P[5])||')'
 
# "!Cycle Statistics"
W_net=W_t_hp+W_t_lp-W_p	#"net work"
Eff=W_net/Q_in	#"cycle eficiency"

print(W_net, Q_out, Eff)
 
