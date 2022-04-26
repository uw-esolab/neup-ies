
# Install CoolProp using pip:
# >> pip install CoolProp

# Author: Mike Wagner | esolab.engr.wisc.edu
# Property lookups are documented at the following URL:
# http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table

import CoolProp.CoolProp as CP

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
        self.wants = "TPHSXD"
        if kwargs:
            if 'T' in kwargs: 
                self.T = kwargs['T']
            if 'P' in kwargs:
                self.P = kwargs['P']
            if 'H' in kwargs:
                self.H = kwargs['H']
            if 'S' in kwargs:
                self.S = kwargs['S']
            if 'X' in kwargs:
                self.X = kwargs['X']
            if 'wants' in kwargs:
                self.wants = kwargs['wants']

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
            if 'H' in self.wants: self.H =  CP.PropsSI("H","T",self.T,"P",self.P, "Water")
            if 'S' in self.wants: self.S =  CP.PropsSI("S","T",self.T,"P",self.P, "Water")
            self.X = -999
        elif self.T and self.H:
            if 'P' in self.wants: self.P = CP.PropsSI("P","H",self.H,"T",self.T, "Water")*1000
            if 'S' in self.wants: self.S =  CP.PropsSI("S","T",self.T,"P",self.P, "Water")
            if 'X' in self.wants: self.X =  CP.PropsSI("Q", "T", self.T, "H", self.H, "Water")
        elif self.P and self.H: 
            if 'T' in self.wants: self.T =  CP.PropsSI("T","H",self.H,"P",self.P, "Water") 
            if 'S' in self.wants: self.S =  CP.PropsSI("S","H",self.H,"P",self.P, "Water")
            if 'X' in self.wants: self.X =  CP.PropsSI("Q", "P", self.P, "H", self.H, "Water")
        elif self.P and hasattr(self, 'X'):
            if 'T' in self.wants: self.T =  CP.PropsSI("T", "P", self.P, "Q", self.X, "Water")
            if 'H' in self.wants: self.H =  CP.PropsSI("H", "P", self.P, "Q", self.X, "Water")
            if 'S' in self.wants: self.S =  CP.PropsSI("S", "P", self.P, "Q", self.X, "Water")
        elif self.T and hasattr(self, 'X'):
            if 'P' in self.wants: self.P = CP.PropsSI("P", "T", self.T, "Q", self.X, "Water")*1000
            if 'H' in self.wants: self.H =  CP.PropsSI("H", "T", self.T, "Q", self.X, "Water")
            if 'S' in self.wants: self.S =  CP.PropsSI("S", "T", self.T, "Q", self.X, "Water")
        elif self.P and self.S:
            if 'T' in self.wants: self.T =  CP.PropsSI("T","P",self.P,"S",self.S, "Water")
            if 'H' in self.wants: self.H =  CP.PropsSI("H","P",self.P,"S",self.S, "Water")
            if 'X' in self.wants: self.X =  CP.PropsSI("Q", "P", self.P, "S", self.S, "Water")
        elif self.T and self.S:
            if 'P' in self.wants: self.P = CP.PropsSI("P","T",self.T,"S",self.S, "Water")*1000
            if 'H' in self.wants: self.H =  CP.PropsSI("H","T",self.T,"S",self.S, "Water")
            if 'x' in self.wants: self.X =  CP.PropsSI("Q", "T", self.T, "S", self.S, "Water")
        else:
            raise RuntimeError("State evaluation failed. Insufficient information provided.")
        
        if 'D' in self.wants: 
            self.rho =  CP.PropsSI("D", "P", self.P, "H", self.H, "Water")
            self.v = 1./self.rho 

def enthalpy(**kwargs): 
    S = State("none", wants='H', **kwargs)
    return S.H

def entropy(**kwargs):
    S = State("none", wants='S', **kwargs)
    return S.S
def volume(**kwargs):
    return State("none", wants='D', **kwargs).v
def density(**kwargs):
    return State("none", wants='D', **kwargs).rho
def temperature(**kwargs):
    return State("none", wants='T', **kwargs).T
def quality(**kwargs):
    return State("none", wants='X', **kwargs).X
def t_sat(P):
    return State("none", P=P, X=0.).T

# ######################################################################################################

