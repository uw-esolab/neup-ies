
# Install CoolProp using pip:
# >> pip install CoolProp

# Author: Mike Wagner | esolab.engr.wisc.edu
# Property lookups are documented at the following URL:
# http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table

import CoolProp.CoolProp as CP
import pint 

un = pint.UnitRegistry()

class State:
    def __init__(self, name : str, **kwargs):
        """
        Parameters
        ----------
        T : (float) 
            [K] temperature
        P : (float)
            [Pa] pressure 
        H : (float)
            [J/kg] enthalpy
        S : (float)
            [J/kg-K] entropy per unit mass 
        X : (float)
            [-] quality
        """
        self.name = name
        self.T = self.P = self.H = self.S = self.v = self.rho = None
        if kwargs:
            if 'T' in kwargs: 
                self.T = un.Quantity(kwargs['T'], 'degK')
            if 'P' in kwargs:
                self.P = un.Quantity(kwargs['P'], 'Pa')
            if 'H' in kwargs:
                self.H = un.Quantity(kwargs['H'], 'J/kg')
            if 'S' in kwargs:
                self.S = un.Quantity(kwargs['S'], 'J/(kg*K)')
            if 'X' in kwargs:
                self.X = kwargs['X']*un.parse_units('')

            self.evaluate()
    def __str__(self):
        xx = -999
        if hasattr(self, 'X'):
            if self.X != None:
                xx = self.X 
        return "{:>20s}> T: {:.1f}\tP: {:.1f}\tH: {:.1f}\tS: {:.1f}\tX: {:.3f}".format(self.name, self.T, self.P, self.H, self.S, xx)

    def evaluate(self): #, method='PH'):
        """
        
        """

        



        if self.T and self.P: 
            self.H = un.Quantity( CP.PropsSI("H","T",self.T.m,"P",self.P.m, "Water"), 'J/kg' )
            self.S = un.Quantity( CP.PropsSI("S","T",self.T.m,"P",self.P.m, "Water"), 'J/(kg*K)' )
            self.X = un.Quantity( -999, '')
        elif self.T and self.H:
            self.P = CP.PropsSI("P","H",self.H.m,"T",self.T.m, "Water")*un.Pa
            self.S = un.Quantity( CP.PropsSI("S","T",self.T.m,"P",self.P.m, "Water"), 'J/(kg*K)' )
            self.X = un.Quantity( CP.PropsSI("Q", "T", self.T.m, "H", self.H.m, "Water"), '')
        elif self.P and self.H: 
            self.T = un.Quantity( CP.PropsSI("T","H",self.H.m,"P",self.P.m, "Water") , 'degK' )
            self.S = un.Quantity( CP.PropsSI("S","H",self.H.m,"P",self.P.m, "Water"), 'J/(kg*K)' )
            self.X = un.Quantity( CP.PropsSI("Q", "P", self.P.m, "H", self.H.m, "Water"), '')
        elif self.P and hasattr(self, 'X'):
            self.T = un.Quantity( CP.PropsSI("T", "P", self.P.m, "Q", self.X.m, "Water"), 'degK' )
            self.H = un.Quantity( CP.PropsSI("H", "P", self.P.m, "Q", self.X.m, "Water"), 'J/kg' )
            self.S = un.Quantity( CP.PropsSI("S", "P", self.P.m, "Q", self.X.m, "Water"), 'J/(kg*K)' )
        elif self.T and hasattr(self, 'X'):
            self.P = CP.PropsSI("P", "T", self.T.m, "Q", self.X.m, "Water")*un.Pa
            self.H = un.Quantity( CP.PropsSI("H", "T", self.T.m, "Q", self.X.m, "Water"), 'J/kg' )
            self.S = un.Quantity( CP.PropsSI("S", "T", self.T.m, "Q", self.X.m, "Water"), 'J/(kg*K)' )
        elif self.P and self.S:
            self.T = un.Quantity( CP.PropsSI("T","P",self.P.m,"S",self.S.m, "Water"), 'degK' )
            self.H = un.Quantity( CP.PropsSI("H","P",self.P.m,"S",self.S.m, "Water"), 'J/kg' )
            self.X = un.Quantity( CP.PropsSI("Q", "P", self.P.m, "S", self.S.m, "Water"), '')
        elif self.T and self.S:
            self.P = CP.PropsSI("P","T",self.T.m,"S",self.S.m, "Water")*un.Pa
            self.H = un.Quantity( CP.PropsSI("H","T",self.T.m,"S",self.S.m, "Water"), 'J/kg' )
            self.X = un.Quantity( CP.PropsSI("Q", "T", self.T.m, "S", self.S.m, "Water"), '')
        else:
            raise RuntimeError("State evaluation failed. Insufficient information provided.")
        
        self.rho = un.Quantity( CP.PropsSI("D", "P", self.P.m, "H", self.H.m, "Water"), 'kg/m^3' )
        self.v = 1./self.rho 

def enthalpy(**kwargs): 
    S = State("none", **kwargs)
    return S.H.m

def entropy(**kwargs):
    S = State("none", **kwargs)
    return S.S.m
def volume(**kwargs):
    return State("none", **kwargs).v.m
def density(**kwargs):
    return State("none", **kwargs).rho.m
def temperature(**kwargs):
    return State("none", **kwargs).T.m
def quality(**kwargs):
    return State("none", **kwargs).X.m
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
 
