from sympy import N
from state import *
from math import sqrt
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
def fwh_TTD( m_dot_D, h_D_in, P_D, m_dot_FW, h_FW_in, P_FW, h_T_in, P_T, DT_TTD, DT, N_hxrs): 
    
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

        
    iter = 1
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
        DT_TTD_c = T_FW_out - (t_sat( P = P_T)+0.01)
        diff_c = abs(DT_TTD - DT_TTD_c)

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
        DT_TTD_d = T_FW_out - (t_sat( P = P_T)+0.01)
        diff_d = abs(DT_TTD - DT_TTD_d)
        #golden section search iteration
        if diff_c < diff_d:
            m_b = m_d
        else:
            m_a = m_c
        m_c = m_b - (m_b -m_a)/gr
        m_d = m_a + (m_b-m_a)/gr
        iter = iter + 1
        #print("iter:" + str(iter) + " -- error: " + str(err))

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


