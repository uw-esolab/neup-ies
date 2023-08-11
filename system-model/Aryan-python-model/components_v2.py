from sympy import N
from state import *
from math import *
import numpy as np

##  Closed Feed water heater ==================================================================

def fwh(m_dot_D,h_D_in,P_D,m_dot_FW,h_FW_in,P_FW,m_dot_T,h_T_in,P_T,DT,N_hxrs):
    """
    Closed feed water heater function. Returns hot and cold temperature arrays.

    Arguments:
        m_dot_D --> Drain mass flow rate [kg/s]
        h_D_in --> Drain inlet enthalpy [J/kg]
        P_D --> Drain flow pressure [Pa]
        m_dot_FW --> Feed water heater mass flow [kg/s]
        h_FW_in --> Feed water heater inlet enthalpy [J/kg]
        P_FW --> Feed water heater pressure [Pa]
        m_dot_T --> Turbine mass flow rate [kg/s]
        h_T_in -->  Turbine flow inlet enthalpy [J/kg]
        P_T --> Turbine flow pressure [Pa]
        N_hxrs --> number of sub heat exchangers
    
    Outputs:
        UA_ds --> desuperheating region conductance [W/K]
        UA_cond --> condensing region conductance [W/K]
        UA_dc --> drain cooler conductance [W/K]
        NTU --> NTU array
        T_H --> hot side temperature array [K]
        T_C --> cold side temperature array [K]        
    """
    # get boundary states
    T_T_in = temperature(P=P_T,H=h_T_in)                            # turbine flow inlet temperature
    T_T_out = t_sat(P=P_T) - 0.01                                   # turbine flow outlet temperature
    h_T_out = enthalpy(P=P_T,T=T_T_out)                             # turbine flow outlet enthalpy
    m_dot_D_int = m_dot_D + m_dot_T                                 # drain internal mass flow rate
    h_D_int = (m_dot_T*h_T_out + m_dot_D*h_D_in)/m_dot_D_int        # drain internal inlet enthalpy
    P_D_int = P_T                                                   # TODO: ask about PD_int. shouldn't it be min(P_D,P_T)
    T_D_int = temperature(P=P_D_int,H=h_D_int)                      # internal drain inlet temperature
    T_FW_in = temperature(H=h_FW_in,P=P_FW)                         # feed water flow inlet temperature
    T_D_out = T_FW_in + DT                                          # drain outlet temperature
    h_D_out = enthalpy(T=T_D_out,P=P_D_int)                         # drain outlet enthalpy
    q_dot_hx = m_dot_D_int*(h_D_int-h_D_out)                        # net heat transfer
    h_FW_int = h_FW_in + q_dot_hx/m_dot_FW                          # feed water flow internal enthalpy
    T_FW_int = temperature(H=h_FW_in,P=P_D_int)                     # feed water flow internal temperature
    q_dot_ds_cond = m_dot_T*(h_T_in - h_T_out)
    h_FW_out = h_FW_int + q_dot_ds_cond/m_dot_FW                    # feed water flow outlet enthalpy
    T_FW_out = temperature(H=h_FW_out,P=P_FW)                       # feed water flow outlet temperature
    q_dot_tot = m_dot_FW*(h_FW_out - h_FW_in)                       # total heat transfer 
    h_T_vap = enthalpy(P=P_T,T=t_sat(P=P_T)+0.01)                   # turbine enthalpy at saturated vapor

    if(h_T_in-h_T_vap > 0):
        q_dot_ds = m_dot_T*(h_T_in - h_T_vap)
        q_dot_cond = m_dot_T*(h_T_vap - h_T_out)
    else:
        q_dot_ds = 0
        q_dot_cond = m_dot_T*(h_T_in - h_T_out)
    
    # build temperature arrays
    nn = 3*N_hxrs+2
    h_H = np.zeros(nn)                                              # hot side enthalpy array
    T_H = np.zeros(nn)                                              # hot side temperature array
    h_C = np.zeros(nn)                                              # cold side enthalpy array
    T_C = np.zeros(nn)                                              # cold side temperature array
    C_dot_H = np.zeros(nn)                                          # hot side capacitance rates
    C_dot_C = np.zeros(nn)                                          # cold side capacitance rates
    q_dot = np.zeros(nn)                                            # heat transfer
    eff = np.zeros(nn)                                              # sub hx effectiveness
    NTU = np.zeros(nn)                                              # sub hx NTU
    UA = np.zeros(nn)                                               # sub hx UA
    DELTA_T = np.zeros(nn)                                          # sub hx approach temperature difference

    h_H[1] = h_T_in
    T_H[1] = T_T_in
    h_C[1] = h_FW_out 
    T_C[1] = T_FW_out 
    # desuperheater
    for i in range(2,(N_hxrs+1)+1) : 	                            # building hot enthalpy and temperature array
        h_H[i] = h_H[i-1]-q_dot_ds/(N_hxrs*m_dot_T)
        T_H[i] = temperature(H = h_H[i],P=P_T)
    for i in range(2,(N_hxrs+1)+1) : 	                            # building cold enthalpy and temperature array
        h_C[i] = h_C[i-1]-q_dot_ds/(N_hxrs*m_dot_FW)	
        T_C[i]=temperature(H = h_C[i],P=P_FW)	
    # condensing
    for i in range((N_hxrs+2),(2*N_hxrs)+1) : 	                    # building hot enthalpy and temperature array
        h_H[i] = h_H[i-1]-q_dot_cond/(N_hxrs*m_dot_T)
        T_H[i] = temperature(H = h_H[i],P=P_T)
    for i in range((N_hxrs+2),(2*N_hxrs+1)+1) : 	                # building cold enthalpy and temperature array
        h_C[i] = h_C[i-1]-q_dot_cond/(N_hxrs*m_dot_FW)	
        T_C[i]=temperature(H = h_C[i],P=P_FW)	
    h_H[2*N_hxrs+1] = h_D_int	                                    # drain addition enthalpy
    T_H[2*N_hxrs+1] =  T_D_int	                                    # drain addition temperature
    # hx
    for i in range((2*N_hxrs+2),(3*N_hxrs+1)+1) : 	                # building hot enthalpy and temperature array
        h_H[i] = h_H[i-1]-q_dot_hx/(N_hxrs*m_dot_D_int)
        T_H[i] = temperature(H = h_H[i],P=P_D_int)
    for i in range((2*N_hxrs+2),(3*N_hxrs+1)+1) : 	                # building cold enthalpy and temperature array
        h_C[i] = h_C[i-1]-q_dot_hx/(N_hxrs*m_dot_FW)	
        T_C[i]=temperature(H = h_C[i],P=P_FW)	

    # Capacitance rates array
    for i in range(1,(3*N_hxrs+1)):
        dTH = (T_H[i] - T_H[i+1])
        dTC = (T_C[i] - T_C[i+1])

        if(i <= 2*N_hxrs):
            if(dTH <= 0):
                C_dot_H[i] = 1e99
            else:
                C_dot_H[i] = m_dot_T * (h_H[i] - h_H[i+1])/dTH
            if(dTC <= 0):
                C_dot_C[i] = 1e99
            else:
                C_dot_C[i] = m_dot_FW*(h_C[i] - h_C[i+1])/dTC
        else:
            if(dTH <= 0):
                C_dot_H[i] = 1e99
            else:
                C_dot_H[i] = (m_dot_T + m_dot_D) * (h_H[i] - h_H[i+1])/dTH
            if(dTC <= 0):
                C_dot_C[i] = 1e99
            else:
                C_dot_C[i] = m_dot_FW*(h_C[i] - h_C[i+1])/dTC

    # Get effectiveness and NTU array
    for i in range(1,(3*N_hxrs+1)):
        q_dot[i] = m_dot_FW*(h_C[i]-h_C[i+1])
        # effectiveness of sub heat exchanger
        cdotdT=min(C_dot_H[i],C_dot_C[i])*(T_H[i]-T_C[i+1])
        # if cdotdT==0.: cdotdT=1e-6
        eff[i] = q_dot[i]/cdotdT

        if(eff[i] <= 0):
            NTU[i] = 0
            UA[i] = 0
        else:
            NTU[i] = NTU_counterflow(eff[i],C_dot_H[i],C_dot_C[i])      # NTU in sub heat exchanger                                          
            UA[i] = NTU[i] *min(C_dot_H[i],C_dot_C[i])                  # conductance in sub heat exchanger
    
    # get net conductance values
    UA_ds = UA[1:N_hxrs+1].sum()  # desuperheater conductance
    UA_cond = UA[N_hxrs+1:2*N_hxrs+1].sum() # condensor conductance
    UA_dc = UA[2*N_hxrs+1:3*N_hxrs+1].sum() # drain cooler conductance
    UA_tot = UA_ds + UA_cond + UA_dc    # total conductance

    return [
        UA_ds, 
        UA_cond, 
        UA_dc,
        NTU[1:3*N_hxrs+2],
        T_H[1:3*N_hxrs+2],
        T_C[1:3*N_hxrs+2], 
    ]

## subcounter flow heat exchanger ==============================================================
# TODO: add subcounter flow hx for on design from Brian's model.
# TODO: ask Mike about Molten Salt model tables/incompressible model for charging and discharging
def sub_counterflow(F_H,m_dot_H,h_H_in,P_H_in,F_C,m_dot_C,h_C_in,P_C_in,DT,N_hxrs):
    """
    Counterflow heat exchanger model for the CSP HX. Currently valid for OFF design.

    Arguments:
        F_H --> hot side fluid (OFF - Water [W], C - Water [W], D - Salt [S])
        m_dot_H --> hot side mass flow rate [kg/s]
        h_H_in --> hot side inlet enthalpy [J/kg]
        P_H_in --> hot side pressure [Pa]
        F_C --> cold side fluid (OFF - Water [W], C - Salt [S], D - Water [W])
        m_dot_C --> cold side mass flow rate [kg/s]
        h_C_in --> cold side inlet enthalpy [J/kg]
        P_C_in --> cold side pressure [Pa]
        DT --> hx approach temperature difference
        N_hxrs --> number of sub heat exchangers
    
    Outputs:
        h_H_out --> hot side outlet enthalpy [J/kg]
        h_C_out --> cold side outlet enthalpy [J/kg]
        q_dot_tot --> net heat transfer [W]
        UA_tot --> net conductance
        T_H --> hot side temperature array [K]
        T_C --> cold side temperature array [K]
    """
    # initialize arrays
    T_H = np.zeros(N_hxrs + 1)
    h_H = np.zeros(N_hxrs + 1)
    T_C = np.zeros(N_hxrs + 1)
    h_C = np.zeros(N_hxrs + 1)
    q_dot = np.zeros(N_hxrs + 1)
    C_dot_C = np.zeros(N_hxrs + 1)
    C_dot_H = np.zeros(N_hxrs + 1)
    eff = np.zeros(N_hxrs + 1)
    NTU = np.zeros(N_hxrs + 1)
    UA = np.zeros(N_hxrs + 1)

    # check fluid type for inlet temperature
    # hot side
    if(F_H == 'W' or F_H.lower() == "water"):
        T_H[1] = temperature(H=h_H_in,P=P_H_in)
    elif(F_H == 'S' or F_H.lower() == "salt"):
        # TODO: add properties for molten salt
        T_H[1] = 0  # place holder
    else:
        raise Exception("Fluid type is not accepted.")
    # cold side
    if(F_C == 'W' or F_C.lower() == "water"):
        T_C[N_hxrs+1] = temperature(H=h_C_in,P=P_C_in)
        h_H_out_min = enthalpy(T=T_C[N_hxrs+1],P=P_C_in)        

    elif(F_C == 'S' or F_C.lower() == "salt"):
        # TODO: add properties for molten salt
        T_C[N_hxrs+1] = 0  # place holder
    else:
        raise Exception("Fluid type is not accepted.")
    
    T_H[N_hxrs+1] = T_C[N_hxrs+1] + DT                                          # assuming approach temperature on cold side
    # hot side outlet enthalpy
    if(F_H == 'W' or F_H.lower() == "water"):
        h_H_out = enthalpy(T=T_H[N_hxrs+1],P=P_H_in)            
    elif(F_H == 'S' or F_H.lower() == "salt"):
        # TODO: add properties for molten salt
        h_H_out = 0  # place holder
    else:
        raise Exception("Fluid type is not accepted.")
    
    q_dot_tot = m_dot_H*(h_H_in-h_H_out)                                        # total heat transfer in HX
    h_C_out = h_C_in + q_dot_tot/m_dot_C                                        # cold side outlet enthalpy

    # cold side outlet temperature
    if(F_C == 'W' or F_C.lower() == "water"):
        T_C[1] = temperature(H=h_C_out,P=P_C_in)
    elif(F_C == 'S' or F_C.lower() == "salt"):
        # TODO: add properties for molten salt
        T_C[1] = 0  # place holder
    else:
        raise Exception("Fluid type is not accepted.")

    h_H[1] = h_H_in                                                             # hot side inlet enthalpy
    h_C[1] = h_C_out                                                            # cold side outlet enthalpy

    for i in range(2,(N_hxrs+2)):
        # hot side
        h_H[i] = h_H[i-1] - q_dot_tot/(N_hxrs*m_dot_H)                          # hot side enthalpy for sub heat exchangers
        # get sub heat exchanger temperatures
        if(F_H == 'W' or F_H.lower() == "water"):
            T_H[i] = temperature(H=h_H[i],P=P_H_in)
        elif(F_H == 'S' or F_H.lower() == "salt"):
            # TODO: add properties for molten salt
            T_H[i] = 0  # place holder
        else:
            raise Exception("Fluid type is not accepted.")
        # cold side
        h_C[i] = h_C[i-1] - q_dot_tot/(N_hxrs*m_dot_C)                          # cold side enthalpy for sub heat exchangers
        # get sub heat exchanger temperatures
        if(F_C == 'W' or F_C.lower() == "water"):
            T_C[i] = temperature(H=h_C[i],P=P_C_in)
        elif(F_C == 'S' or F_C.lower() == "salt"):
            # TODO: add properties for molten salt
            T_C[i] = 0  # place holder
        else:
            raise Exception("Fluid type is not accepted.")
        
        q_dot_max = (h_H[1] - h_H_out_min)*m_dot_H

    # Checking hot end pinch point
    if(T_H[1] - T_C[1] < DT):
        T_C[1] = T_H[1] - DT                                                    # switch pinch point to hot end 
        # cold side outlet temperature
        if(F_C == 'W' or F_C.lower() == "water"):
            h_C_out = enthalpy(T=T_C[1],P=P_C_in)
        elif(F_C == 'S' or F_C.lower() == "salt"):
            # TODO: add properties for molten salt
            h_C_out = 0  # place holder
        else:
            raise Exception("Fluid type is not accepted.")
        q_dot_tot = m_dot_C*(h_C_out-h_C_in)                                    # total heat transfer in HX
        h_H_out = h_H_in - q_dot_tot/m_dot_H                                    # hot side outlet enthalpy

        # hot side outlet temperature
        if(F_H == 'W' or F_H.lower() == "water"):
            T_H[N_hxrs+1] = temperature(H=h_H_out,P=P_H_in)
            h_C_out_max = enthalpy(T=T_H[1],P=P_C_in)
        elif(F_H == 'S' or F_H.lower() == "salt"):
            # TODO: add properties for molten salt
            T_H[N_hxrs+1] = 0  # place holder
        else:
            raise Exception("Fluid type is not accepted.")
        
        h_H[1] = h_H_in                                                         # hot side inlet enthalpy
        h_C[1] = h_C_out                                                        # cold side outlet enthalpy

        for i in range(2,(N_hxrs+2)):
            # hot side
            h_H[i] = h_H[i-1]-q_dot_tot/(N_hxrs*m_dot_H)                        # hot side enthalpy for sub hx
            if(F_H == 'W' or F_H.lower() == "water"):
                T_H[i] = temperature(H=h_H[i],P=P_H_in)
            elif(F_H == 'S' or F_H.lower() == "salt"):
                # TODO: add properties for molten salt
                T_H[i] = 0  # place holder
            else:
                raise Exception("Fluid type is not accepted.")
            # cold side
            h_C[i] = h_C[i-1]-q_dot_tot/(N_hxrs*m_dot_C)                        # cold side enthalpy for sub hx
            if(F_C == 'W' or F_C.lower() == "water"):
                T_C[i] = temperature(H=h_C[i],P=P_C_in)
            elif(F_C == 'S' or F_C.lower() == "salt"):
                # TODO: add properties for molten salt
                T_C[i] = 0  # place holder
            else:
                raise Exception("Fluid type is not accepted.")
        q_dot_max = m_dot_C*(h_C_out_max - h_C[N_hxrs+1]) 

    for i in range(1,(N_hxrs+1)):
        # get capacitance rates
        C_dot_H[i] = m_dot_H*(h_H[i]-h_H[i+1])/(T_H[i]-T_H[i+1])
        C_dot_C[i] = m_dot_C*(h_C[i]-h_C[i+1])/(T_C[i]-T_C[i+1])

        q_dot[i] = m_dot_C*(h_C[i]-h_C[i+1])                                    # heat transfer of sub hx
        eff[i] = q_dot[i]/(min(C_dot_H[i],C_dot_C[i])*(T_H[i]-T_C[i+1]))        # effectiveness of sub hx
        NTU[i] = NTU_counterflow(eff[i],C_dot_H[i],C_dot_C[i])                  # NTU required by sub hx 
        UA[i] = NTU[i]*min(C_dot_H[i],C_dot_C[i])                               # conductance in sub hx

    UA_tot = sum(UA)                                                            # total conductance     

    return [
        h_H_out,
        h_C_out,
        q_dot_tot,
        UA_tot,
        T_H[1:(N_hxrs+1)],
        T_C[1:(N_hxrs+1)],
    ]

## counter flow heat exchanger procedure (from EES) ============================================
def NTU_counterflow(eff,C_dot_H,C_dot_C):
    """
    Counter flow hx model. Returns NTU based on effectiveness.

    Arguments:
        eff --> effectiveness of heat exchanger
        C_dot_H --> hot side capacitance rate
        C_dot_C --> cold side capacitance rate

    Output:
        NTU --> number of tranfer units (Eqn from Table 8-2 -- Heat Transfer, Klein and Nellis)    
    """
    C_dot_max = max(C_dot_H,C_dot_C)
    C_dot_min = min(C_dot_H,C_dot_C)
    CR = C_dot_min/C_dot_max

    if(CR>0.99):
        NTU = eff/(1-eff)
    else:
        NTU = log((1-eff*CR)/(1-eff))
    # else:
    #     raise Exception("Capacitance Ratio is greater than 1.")

    return NTU

##  Open Feed Water Heater =====================================================================
def open_fwh( m_dot_fw, h_fw, P_fw, m_dot_turb, h_turb, P_turb, m_dot_drain, h_drain, P_drain):
    """Open feed water heater model. Solves for pressure and enthalpy out."""

    #F$ = working fluid
    #fw = feed water | turb = turbine | drain = drain | out = outlet
    #m_dot_ = mass flow [kg/s]
    #h_ = enthalpy [J/kg]
    #P_ = pressure [Pa]

    m_dot_out = m_dot_fw + m_dot_turb + m_dot_drain
    P_out = (P_turb*m_dot_turb+P_fw*m_dot_fw+P_drain*m_dot_drain)/m_dot_out
    h_out = (m_dot_fw*h_fw + m_dot_turb*h_turb + m_dot_drain*h_drain)/m_dot_out
    return m_dot_out, h_out, P_out

# =============================================================================================

##  Turbine ===================================================================================

## New turbine procedure that accounts for splitting of mass flow to extraction
def turbine(m_dot,h_in,P_in,P_out,eta,Q_dot_fwh,h_fwh_out,m_dot_D,h_D):
    """Turbine Procedure that accounts for extraction with splitting and mixing."""   

    s_in= entropy(P=P_in,H=h_in)                # Inlet entropy
    h_out_s= enthalpy(P=P_out,S=s_in)           # isentropic exit enthalpy
    h_out=h_in-eta*(h_in-h_out_s)               # actual outlet enthalpy
    W_dot=m_dot*(h_in-h_out)                    # power   
    x_out = quality(P=P_out,H=h_out)            # outlet quality 

    if(x_out<1 and x_out>0):
        # splitter
        m_v_in = m_dot * x_out                  # vapor into splitter
        m_l_in = m_dot * (1-x_out)              # liquid into splitter
        m_l_out = m_l_in * 0.7                  # liquid pulled out of flow 
        x_l_out = 0
        m_v_out = m_v_in + m_l_in - m_l_out     # mixed mass flow out of splitter
        x_v_out = m_v_in/m_v_out                # quality out of splitter mix flow
        h_v_out = enthalpy(P=P_out,X=x_v_out)   # enthalpy out of splitter mass flow
        h_l_out = enthalpy(P=P_out,X=x_l_out)

        # mixer
        # explicit solution for mixer inlet
        m_m_in = (Q_dot_fwh - m_l_out*h_l_out + m_l_out*h_fwh_out + m_dot_D*h_fwh_out - m_dot_D*h_D)/(h_v_out - h_fwh_out)
        m_fwh = m_l_out + m_m_in                # mass flow of turbine exhaust
        m_t = m_dot - m_fwh                     # mass flow to next turbine
        h_t_in = h_v_out                        # enthalpy to next turbine
        # enthalpy outlet fwh
        h_fwh_in = (Q_dot_fwh + (m_fwh + m_dot_D)*h_fwh_out -m_dot_D*h_D)/m_fwh
        x_t_in = quality(H=h_t_in,P=P_out)      # quality to next turbine
    else:
        h_fwh_in = h_out                        # enthalpy outlet to fwh
        h_t_in = h_out                          # enthalpy to next turbine
        x_t_in = x_out                          # quality to next turbine
        m_fwh = (Q_dot_fwh - m_dot_D*h_D + h_fwh_out*m_dot_D)/(h_fwh_in - h_fwh_out)
        m_t = m_dot - m_fwh

    return h_fwh_in, h_t_in, m_fwh, m_t, x_t_in, W_dot


## simple turbine with no split
def turbine_simple(m_dot, h_in, P_in, P_out, eta):
    """Simple Turbine model with isentropic efficiency. Assumes no mixing."""

    #F$ = working fluid
    #m_dot = mass flow rate [kg/s]
    #h_in = inlet enthalpy [J/kg]
    #P_in = inlet pressure [Pa]
    #P_out = outlet pressure [Pa]
    #eta = isentropic efficiency [-]
    #h_out = outlet enthalpy [J/kg]
    #W_out = turbine work [W]

    s_in=entropy(P=P_in,H=h_in)    #Inlet entropy
    h_out_s=enthalpy(P=P_out,S=s_in)    #isentropic exit enthalpy
    h_out=h_in-eta*(h_in-h_out_s)    #actual outlet enthalpy
    W_dot=m_dot*(h_in-h_out)    #power    

    return h_out, W_dot
 

##  Condenser ===================================================================================
def condenser( m_dot, h_in, P_in, T_r, DT):
    """Simple condensor model with supplied approach temperature difference. Returns Q_dot and effectiveness."""

    T_in=temperature(H=h_in,P=P_in)             # inlet temperature for real fluid
    if (T_in>T_r):
        T_out=T_r+DT                            # outlet temperature
    else:
        T_out=T_r-DT                            # outlet temperature
    P_out=P_in                                  # outlet pressure
    h_out=enthalpy(T=T_out,P=P_out)             # outlet enthalpy
    h_out_i=enthalpy(T=T_r,P=P_out)             # ideal outlet enthalpy

    if (T_in>T_r):    
        Q_dot=m_dot*(h_in-h_out)                # heat transfer out of fluid
    else:
        Q_dot=m_dot*(h_out-h_in)                # heat transfer into fluid
    eff=(h_in-h_out)/(h_in-h_out_i)             # effectiveness
 
    return h_out, P_out, Q_dot, eff

#Pumps ==========================================================================================
 
def pump( m_dot, h_in, P_in, P_out, eta):
    """Simple Pump model with isentropic efficiency. Returns oulet enthalpy and required power."""

    v_in = volume( P = P_in, H = h_in)          # inlet specific volume
    DELTA_P = P_out-P_in                        # Pump pressure drop
    W_dot_P = ((m_dot*v_in)*DELTA_P)/eta        # Power supplied to pump
    h_out = h_in + W_dot_P/m_dot                # Pump outlet enthalpy
    return h_out, W_dot_P