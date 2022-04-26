from sympy import N
from state import *
from math import sqrt
import numpy as np

# #Sub Heat Exchanger ====================================================================
 
def fwh( m_dot_D, h_D_in, P_D, m_dot_FW, h_FW_in, P_FW, h_T_in, P_T, DT_int, DT, N_hxrs): 
    
    # m_T, h_FW_out, h_D_out,T_H[1..3*N_hxrs+1], T_C[1..3*N_hxrs+1])
	
    # $VarInfo h_H[] units='J/kg'
	# $VarInfo h_C[] units='J/kg'
	# $VarInfo T_H[] units='K' lower = 273 upper = 5000
	# $VarInfo T_C[] units='K' lower = 273 upper = 5000
	# $VarInfo q_dot[] units='W'
	# $VarInfo DELTA_T[] units='K'
	# $VarInfo m_T units='kg/s' lower = 1 upper = 100
 
	# $VarInfo q_dot_tot units='W' lower = 0 
	# $VarInfo q_dot_ds units='W' lower = 0 
	# $VarInfo q_dot_hx units='W' lower = 0 
	# $VarInfo q_dot_cond units='W' lower = 0 

    # Calculate the minimum enthalpy
    # h_min = enthalpy(T=t_sat(2000), X=0.)
    # # Calculate the maximum mass flow rate that can be condensed while keeping above the min. enthalpy
    # # h_T_vap = enthalpy( P = P_T, T = (t_sat( P = P_T)+0.01))
    # h_T_vap = enthalpy( P = P_T, X = 1)
    # h_T_out = enthalpy( T = T_T_out, P = P_T)
    # q_dot_cond = m_c*(h_T_vap-h_T_out)
    # h_H[1] = h_T_in
    # T_H[1] = T_T_in
    # h_H[i] = h_H[i-1]-q_dot_cond/(N_hxrs*m_c)

        

    m_a = .1 #[kg/s]	#lower limit on mass flow
    m_b = 100 #[kg/s]	#upper limit on mass flow
    tol = 0.0001 #[kg/s]	#convergence tolerance
    err = tol*999
    gr = (1+sqrt(5))/2	#golden ratio
    m_c = m_b - (m_b-m_a)/gr	#middle sections for golden section search
    m_d = m_a + (m_b-m_a)/gr
    while(err>tol):
        err = abs(m_b-m_a)
        #function c --------------------------- 
        T_T_in = temperature( P = P_T, H = h_T_in)
        T_T_out = t_sat( P = P_T) - 0.01 #[K]
        h_T_out = enthalpy( T = T_T_out, P = P_T)
        m_dot_D_int = m_c+m_dot_D
        h_D_int = (m_c*h_T_out+m_dot_D*h_D_in)/m_dot_D_int

        if P_D > 0:
            P_D_int = min(P_D, P_T)
        else:
            P_D_int = P_T

        T_D_int = temperature( H = h_D_int, P = P_D_int)
        T_FW_in = temperature( H = h_FW_in, P = P_FW)
        T_D_out = T_FW_in+DT
        h_D_out = enthalpy( T = T_D_out, P = P_D_int)
        q_dot_hx = m_dot_D_int*(h_D_int-h_D_out)
        h_FW_int = h_FW_in + q_dot_hx/m_dot_FW
        # T_FW_int = temperature( H = h_FW_int, P = P_FW)
        q_dot_ds__cond = m_c*(h_T_in-h_T_out)
        h_FW_out = h_FW_int+q_dot_ds__cond / m_dot_FW
        T_FW_out = temperature( H = h_FW_out, P = P_FW)
        # q_dot_tot = m_dot_FW*(h_FW_out - h_FW_in)
        
        h_T_vap = enthalpy( P = P_T, T = (t_sat( P = P_T)+0.01))
        h_T_vap = min(h_T_vap, h_T_in) #correct this if entering with quality< 1 <<<< MJW
        q_dot_cond = m_c*(h_T_vap-h_T_out)
        
        if h_T_in - h_T_vap > 0:
            q_dot_ds = m_c*(h_T_in - h_T_vap)
        else:
            q_dot_ds = 0 #[W]
        
        # {q_dot_ds = m_c*(h_T_in - h_T_vap)}
        #ds
        nn = 3*N_hxrs+2
        h_H = np.zeros(nn)
        T_H = np.zeros(nn)
        h_C = np.zeros(nn)
        T_C = np.zeros(nn)
        DELTA_T = np.zeros(nn)

        h_H[1] = h_T_in
        T_H[1] = T_T_in
        h_C[1] = h_FW_out 
        T_C[1] = T_FW_out 
        #ds
        for i in range(2,(N_hxrs+1)+1) : 	#building hot enthalpy and temperature array
            h_H[i] = h_H[i-1]-q_dot_ds/(N_hxrs*m_c)
            T_H[i] = temperature(H = h_H[i],P=P_T)
        for i in range(2,(N_hxrs+1)+1) : 	#building cold enthalpy and temperature array
            h_C[i] = h_C[i-1]-q_dot_ds/(N_hxrs*m_dot_FW)	
            T_C[i]=temperature(H = h_C[i],P=P_FW)	
        #cond
        for i in range((N_hxrs+2),(2*N_hxrs)+1) : 	#building hot enthalpy and temperature array
            h_H[i] = h_H[i-1]-q_dot_cond/(N_hxrs*m_c)
            T_H[i] = temperature(H = h_H[i],P=P_T)
        for i in range((N_hxrs+2),(2*N_hxrs+1)+1) : 	#building cold enthalpy and temperature array
            h_C[i] = h_C[i-1]-q_dot_cond/(N_hxrs*m_dot_FW)	
            T_C[i]=temperature(H = h_C[i],P=P_FW)	
        h_H[2*N_hxrs+1] = h_D_int	#drain addition enthalpy
        T_H[2*N_hxrs+1] =  T_D_int	#drain addition temperature
        #hx
        for i in range((2*N_hxrs+2),(3*N_hxrs+1)+1) : 	#building hot enthalpy and temperature array
            h_H[i] = h_H[i-1]-q_dot_hx/(N_hxrs*m_dot_D_int)
            T_H[i] = temperature(H = h_H[i],P=P_D_int)
        for i in range((2*N_hxrs+2),(3*N_hxrs+1)+1) : 	#building cold enthalpy and temperature array
            h_C[i] = h_C[i-1]-q_dot_hx/(N_hxrs*m_dot_FW)	
            T_C[i]=temperature(H = h_C[i],P=P_FW)	
        for i in range(1,(3*N_hxrs+1)+1) : 
            DELTA_T[i] = T_H[i] - T_C[i]
        pinch_c = min(DELTA_T[1:2*N_hxrs+2+1])
        diff_c = abs(DT_int - pinch_c)
        #function d ---------------------------------------------------------- 
        T_T_in = temperature( P = P_T, H = h_T_in)
        T_T_out = t_sat( P = P_T) - 0.01 #[K]
        h_T_out = enthalpy( T = T_T_out, P = P_T)
        m_dot_D_int = m_d+m_dot_D
        h_D_int = (m_d*h_T_out+m_dot_D*h_D_in)/m_dot_D_int
        if P_D > 0:
            P_D_int = min(P_D, P_T)
        else:
            P_D_int = P_T
        
        T_D_int = temperature( H = h_D_int, P = P_D_int)
        T_FW_in = temperature( H = h_FW_in, P = P_FW)
        T_D_out = T_FW_in+DT
        h_D_out = enthalpy( T = T_D_out, P = P_D_int)
        q_dot_hx = m_dot_D_int*(h_D_int-h_D_out)
        h_FW_int = h_FW_in + q_dot_hx/m_dot_FW
        # T_FW_int = temperature( H = h_FW_int, P = P_FW)
        q_dot_ds__cond = m_d*(h_T_in-h_T_out)
        h_FW_out = h_FW_int+q_dot_ds__cond / m_dot_FW
        T_FW_out = temperature( H = h_FW_out, P = P_FW)
        # q_dot_tot = m_dot_FW*(h_FW_out - h_FW_in)
        
        h_T_vap = enthalpy( P = P_T, T = (t_sat( P = P_T)+0.01))
        h_T_vap = min(h_T_vap, h_T_in) #correct this if entering with quality< 1 <<<< MJW
        q_dot_cond = m_d*(h_T_vap-h_T_out)
        
        if h_T_in - h_T_vap > 0:
            q_dot_ds = m_c*(h_T_in - h_T_vap)
        else:
            q_dot_ds = 0 #[W]
        # {
        # q_dot_ds = m_d*(h_T_in - h_T_vap)}
        #ds
        h_H[1] = h_T_in
        T_H[1] = T_T_in
        h_C[1] = h_FW_out 
        T_C[1] = T_FW_out 
        #ds
        for i in range(2,(N_hxrs+1)+1) : 	#building hot enthalpy and temperature array
            h_H[i] = h_H[i-1]-q_dot_ds/(N_hxrs*m_d)
            T_H[i] = temperature(H = h_H[i],P=P_T)
        for i in range(2,(N_hxrs+1)+1) : 	#building cold enthalpy and temperature array
            h_C[i] = h_C[i-1]-q_dot_ds/(N_hxrs*m_dot_FW)	
            T_C[i]=temperature(H = h_C[i],P=P_FW)	
        #cond
        for i in range((N_hxrs+2),(2*N_hxrs)+1) : 	#building hot enthalpy and temperature array
            h_H[i] = h_H[i-1]-q_dot_cond/(N_hxrs*m_d)
            T_H[i] = temperature(H = h_H[i],P=P_T)
        for i in range((N_hxrs+2),(2*N_hxrs+1)+1) : 	#building cold enthalpy and temperature array
            h_C[i] = h_C[i-1]-q_dot_cond/(N_hxrs*m_dot_FW)	
            T_C[i]=temperature(H = h_C[i],P=P_FW)	
        h_H[2*N_hxrs+1] = h_D_int	#drain addition enthalpy
        T_H[2*N_hxrs+1] =  T_D_int	#drain addition temperature
        #hx
        for i in range((2*N_hxrs+2),(3*N_hxrs+1)+1) : 	#building hot enthalpy and temperature array
            h_H[i] = h_H[i-1]-q_dot_hx/(N_hxrs*m_dot_D_int)
            T_H[i] = temperature(H = h_H[i],P=P_D_int)
        for i in range((2*N_hxrs+2),(3*N_hxrs+1)+1) : 	#building cold enthalpy and temperature array
            h_C[i] = h_C[i-1]-q_dot_hx/(N_hxrs*m_dot_FW)	
            T_C[i]=temperature(H = h_C[i],P=P_FW)	
        for i in range(1,(3*N_hxrs+1)+1): 
            DELTA_T[i] = T_H[i] - T_C[i]
        pinch_d = min(DELTA_T[1:2*N_hxrs+2+1])
        diff_d = abs(DT_int - pinch_d)
        #golden section search iteration
        if diff_c < diff_d:
            m_b = m_d
        else:
            m_a = m_c
        m_c = m_b - (m_b -m_a)/gr
        m_d = m_a + (m_b-m_a)/gr

    m_T = (m_b+m_a)/2
    # print("Pinch: {:.2f}:{:.2f} ".format(pinch_c,pinch_d), end="")
    return [
        m_T, 
        h_FW_out, 
        h_D_out,
        T_H[1:3*N_hxrs+1+1], 
        T_C[1:3*N_hxrs+1+1],
        # pinch_c, 
        # pinch_d,
    ]

# ===============================================================================================


#Turbine =============================================
def turbine( m_dot, h_in, P_in, P_out, eta):
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
 
 
#Open Feed Water Heater ================================
def open_fwh( m_dot_fw, h_fw, P_fw, m_dot_turb, h_turb, P_turb, m_dot_drain, h_drain, P_drain):
    #F$ = working fluid
    #fw = feed water | turb = turbine | drain = drain | out = outlet
    #m_dot_ = mass flow [kg/s]
    #h_ = enthalpy [J/kg]
    #P_ = pressure [Pa]
    m_dot_out = m_dot_fw + m_dot_turb + m_dot_drain
    P_out = (P_turb*m_dot_turb+P_fw*m_dot_fw+P_drain*m_dot_drain)/m_dot_out
    h_out = (m_dot_fw*h_fw + m_dot_turb*h_turb + m_dot_drain*h_drain)/m_dot_out
    return m_dot_out, h_out, P_out
 
#Condenser ===========================================
def condenser( m_dot, h_in, P_in, T_r, DT):
    T_in=temperature(H=h_in,P=P_in)    #Inlet temperature for real
    if (T_in>T_r):
        T_out=T_r+DT    #outlet temperature
    else:
        T_out=T_r-DT    #outlet temperature
    
    P_out=P_in    #outlet pressure

    h_out=enthalpy(T=T_out,P=P_out)    #outlet enthalpy
    h_out_i=enthalpy(T=T_r,P=P_out)    #ideal outlet enthalpy

    if (T_in>T_r):    
        Q_dot=m_dot*(h_in-h_out)    #heat transfer out of fluid
    else:
        Q_dot=m_dot*(h_out-h_in)    #heat transfer into fluid
    eff=(h_in-h_out)/(h_in-h_out_i)    #effectiveness
 
    return h_out, P_out, Q_dot, eff

#Pumps =============================================
 
def pump( m_dot, h_in, P_in, P_out, eta):
    v_in = volume( P = P_in, H = h_in)
    DELTA_P = P_out-P_in
    W_dot_P = ((m_dot*v_in)*DELTA_P)/eta
    h_out = h_in + W_dot_P/m_dot
    return h_out, W_dot_P
 

 























if __name__ == "__main__":
     
    # //constants for testing code
    m_given = [0, 619.385, 619.385,68.802,15.98,670.831,670.831,522.217,39.407,438.914,\
    43.896,438.914,19.4,24.018, 28.702,366.794,25.599,12.993,17.286,13.574,297.343,\
    29813.6,29813.6,395.497,395.497,395.497,395.497,395.497,395.497,522.217,522.217,\
    522.217,522.217,522.217,39.407,83.303,102.702,25.599,38.592,55.878,69.452,28.702]
    T_celc = [0, 25,330.41,25,600,394,127,571,364.23,304.98,304.98,569,459.68,362.68,\
    362.68,362.68,253.53,128.72,88.19,60.94,35.79,25,30.79,35.79,35.94,58.16,85.41,\
    105.06,143.24,180.06,185.59,211.95,253.27,284.05,258.82,217.5,191.14,110.61,90.96,\
    63.71,41.49, 39.68]
    T_given = np.array(T_celc)+273.15

    P_bar = [0, 1.002,1,1,1,0.998,0.978,254,67.97,43,43,41.1,20.58,10.44,10.44,10.44,4.374,\
    1.333,0.655,0.208,0.059,1,1,0.059,17.24,15.84,14.54,12.94,11.44,10.04,308.7,303.7,298.4,\
    293.5,64.97,41,19.08,4.174,1.333,0.655,0.208,0.073]
    
    P_given = np.array(P_bar)*100000.
    h_given = np.zeros(len(P_given))
    for i in range(1,len(h_given)):
        h_given[i] = enthalpy(P=P_given[i], T=T_given[i])
    # F$ = 'water'
    
    
    # //testing reference hx \ example of how to use ====================
    
    # //drain inlet properties:

    # m[34] = m_given[34]
    # h[34] = h_given[34]
    # P[34] = P_given[34]
    # //feedwater inlet properties
    # m[31] = m_given[31]
    # h[31] = h_given[31]
    # P[31] =  P_given[31]
    # //turbine extraction properties
    # {m_dot_T = m[10]}
    # h[10] = h_given[10]
    # P[10] = P_given[10]
    # //approach temperatures
    DT_int_2 = 3.0 #[K]	//internal pinch point
    DT_2 = 5.55 #[K]	//approach temperature on cold side
    N_hxrs = 5
    
    ret = fwh(m_given[34], h_given[34], P_given[34], m_given[31], h_given[31], P_given[31], h_given[10],  P_given[10], DT_int_2, DT_2, N_hxrs)
    #: m[10], h[32], h[35], T_H2[1..3*N_hxrs#+1], T_C2[1..3*N_hxrs#+1])
    print(ret)
