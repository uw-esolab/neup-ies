from state import *
from components_v2 import *
import numpy as N
import time
from math import *
from scipy.optimize import root_scalar

##  OFF Design System Model for Supercritical Rankine Steam Cycle =================================
##  Based off of Brian's EES model --> rankine-on-design ==========================================
# --------------------------

class Cycle:
    def __init__(self, mode='O', logfile=None):
        """
        PARAMETERS:
        ===========
        mode    | Select mode of operation --> O = OFF, C = 100% Charging, D = 100% Discharging
        logfile | If path is specified, the system state will be recorded in a log file during iteration
        """
        self.mode = mode 
        self.logfile = logfile 

        self.N_hxrs = 20 #[-]                       # Number of sub heat exchangers
        self.ns = 33                                # number of states used + 1

        # Initialize members
        self.m = N.zeros(self.ns) * float('NaN')    # mass flow array
        self.h = N.zeros(self.ns) * float('NaN')    # enthalpy array
        self.P = N.zeros(self.ns) * float('NaN')    # pressure array
        self.T = N.zeros(self.ns) * float('NaN')    # temperature array
        self.x = N.zeros(self.ns) * float('NaN')    # quality array

    def __str__(self):
        return \
            f'Power (H/I/L): {self.W_dot_gen:.1f} MW ({self.W_dot_HPT:.1f}/{self.W_dot_IPT:.1f}/{self.W_dot_LPT})\n'+\
            f'Efficiency:    {self.eta_cycle:.4f}\n'+\
            f'UA (A/B/C/D/E/F): {self.UA_table["A"][-1]:.0f}/{self.UA_table["B"][-1]:.0f}/{self.UA_table["C"][-1]:.0f}/{self.UA_table["D"][-1]:.0f}/{self.UA_table["E"][-1]:.0f}/{self.UA_table["F"][-1]:.0f}'

    # Define function for logging and reporting data
    # --------------------------
    def __log_header(self,N,names):
        """Writes the headers for the state output file."""
        hh = []
        for name in names:
            h = []
            for i in range(1,N):
                h.append("{:s}{:d}".format(name,i))
            hh.append(h)
        return ','.join([','.join(h) for h in hh]) + '\n'

    def __log_state(self,*args):
        """Returns the state properties for the converged system."""
        # T,P,s,h,x,m
        return ','.join([','.join([str(v) for v in arr[1:]]) for arr in args]) + '\n'

    def solve_cycle(self):
        tstart = time.time()

        # Model start

        #-----------------------------------------------#
        # TODO: update values
        DT = 5.55 #[K]                                                      # Target Approach Temperature 
        #-----------------------------------------------#
        # Start writing to iteration log file 
        if self.logfile:
            fout = open('iter-log2.csv','w')
            fout.write(self.__log_header(self.ns, "T P h x m".split(' ')))                  # provide headers to log file

        # m_s = N.zeros(ns)                                                   # salt mass flow array

        # T_H_A = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX A
        # T_C_A = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX A
        # T_H_B = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX B
        # T_C_B = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX B
        # T_H_C = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX C
        # T_C_C = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX C
        # T_H_D = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX D
        # T_C_D = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX D
        # T_H_E = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX E
        # T_C_E = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX E
        # T_H_F = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX F
        # T_C_F = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX F
        # T_H_G = N.zeros(3*self.N_hxrs+2)                                         # hot side temperature at HX G
        # T_C_G = N.zeros(3*self.N_hxrs+2)                                         # cold side temperature at HX G

        loc = N.zeros(self.N_hxrs*3+2)                                           # position array across each hx
        # Populate location array based on number of sub heat exchangers
        for i in range(1,(3*self.N_hxrs+1)+1):
            loc[i] = (i-1)/(3*self.N_hxrs)

        # Initialize heat exchanger parameters
        if (self.mode.lower() == "o" or self.mode.lower() == "off"):
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
        elif (self.mode.lower() == "c" or self.mode.lower() == "charging"):
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
            # m_s[7] = 190 #[kg/s]                                            # salt mass flow through W2S
        elif (self.mode.lower == "d" or self.mode.lower() == "discharging"):
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
            # m_s[4] = 525 #[kg/s]                                            # salt mass flow through S2W
            self.m[33] = 100 #[kg/s]                                             # steam mass flow through W2S
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
        self.P[14] = 5000 #[Pa]                                                 # Inlet pressure of LP Turbine
        self.P[16] = 2.544e6 #[Pa]                                               # Outlet pressure of LP Turbine

        # Get Mass Flow Rates ==========
        
        # LFR Boiler
        self.P[2] = 3.3e7                                                        # Boiler outlet pressure
        self.T[2] = 632+273.15                                                   # Boiler outlet temperature
        self.h[2] = enthalpy( T = self.T[2], P = self.P[2])                                # Boiler outlet enthalpy
        self.P[1] = self.P[2]                                                         # Boiler inlet pressure, assuming no pressure loss across boiler
        self.T[1] = 340+273.15                                                   # Boiler inlet temperature --> KNOWN
        self.h[1] = enthalpy( T = self.T[1], P = self.P[1])                                # Boiler inlet enthalpy
        Q_dot_LFR = 950e6                                                   # Heat transfer from LFR

        self.m[1] = Q_dot_LFR/(self.h[2]-self.h[1])                                        # Get cycle mass flow rate
        self.m[2] = self.m[1]
        self.h[32] = self.h[31]
        self.m[21] = self.m[1]

        # High Pressure Turbines

        # HPT 1
        self.P[3] = 2.365e7 #[Pa]                                                # 

        # Initial guess for self.m[15]
        self.m[15] = 300 #[kg/s]

        it_max = 150
        tol = 1e-6

        # while((abs(err_iter) > tol and it<it_max) or it<3 ):
        def cycle_eval(m15):
            self.m[15] = m15
            
            # print("Iter {:d}: ".format(it), end='')

            # Temperature Estimation out of LP Pump
            self.P[15] = self.P[14]                                                       # No Pressure Drop across Condenser
            self.T[15] = temperature(P=self.P[15],X=0)                                    # LP Pump inlet temperature
            self.h[15] = enthalpy(X=0, P=self.P[15])                                   # LP Pump inlet enthalpy
            self.h[16], W_dot_LP = pump( self.m[15], self.h[15], self.P[15], self.P[16], eta_LP)
            self.T[16] = temperature(P=self.P[16],H=self.h[16])                                # LP Pump outlet temperature
            self.m[16] = self.m[15]
            # Get FWH temperature rises assuming even temp rise across FWH
            n_fwh = 8                                                           # Number of feed water heaters
            T_rise = (self.T[1] - self.T[16])/n_fwh
            self.T[17] = self.T[16] + T_rise
            self.T[18] = self.T[17] + T_rise
            self.T[19] = self.T[18] + T_rise
            self.T[20] = self.T[19] + T_rise
            self.T[22] = self.T[20] + T_rise
            self.T[23] = self.T[22] + T_rise
            self.T[24] = self.T[23] + T_rise
            self.T[32] = self.T[1]
            # Pressure in FW
            self.P[17] = self.P[16]                                                       # Low Pressure FW
            self.P[18] = self.P[17]
            self.P[19] = self.P[18]
            self.P[20] = self.P[19]

            self.P[22] = self.P[1]                                                        # High Pressure FW
            self.P[23] = self.P[22]
            self.P[24] = self.P[23]
            self.P[32] = self.P[24]

            # Heat Transfer required in hp FW
            for i in range(18,25):
                if i != 21:
                    self.h[i] = enthalpy(T=self.T[i],P=self.P[i])
            self.h[32] = enthalpy(P=self.P[32],T=self.T[32])

            # HP Turbines
            self.m[17] = self.m[16]
            self.m[18] = self.m[17]
            self.m[19] = self.m[18]
            self.m[20] = self.m[19]
            self.m[22] = self.m[21]
            self.m[23] = self.m[22]
            self.m[24] = self.m[23]

            Q_dot_HXE = self.m[22]*(self.h[23]-self.h[22])
            Q_dot_HXF = self.m[23]*(self.h[24]-self.h[23])
            Q_dot_HXG = self.m[24]*(self.h[32]-self.h[24])

            self.P[6] = pressure(T = (self.T[24]+TTD_F), X = 0.1)
            self.P[5] = (self.P[6]+self.P[3])/2
            self.P[7] = pressure(T = (self.T[23]+TTD_E), X = 0.1)
            self.P[8] = pressure(T = self.T[22], X = 0.1)
            self.P[9] = self.P[8]
            self.P[10] = pressure(T = (self.T[20]+TTD_D), X = 0.1)
            self.P[11] = pressure(T = (self.T[19]+TTD_C), X = 0.1)
            self.P[12] = pressure(T = (self.T[18]+TTD_B), X = 0.1)
            self.P[13] = pressure(T = (self.T[17]+TTD_A), X = 0.1)

            # HP energy balance for turbine extraction mass flows

            # HPT 1
            # Define HPT1 states for OFF design
            self.P[4] = self.P[3]
            self.P[25] = self.P[4]
            self.T[25] = self.T[24]+DT_HX_G
            self.h[25] = enthalpy(T = self.T[25], P = self.P[25])
            self.h[4], h_HPT1, self.m[4], m_HPT1, x_HPT1, W_dot_HPT1 = turbine(self.m[2], self.h[2], self.P[2], self.P[3], eta_HPT, Q_dot_HXG, self.h[25], 0, 0)
            Q_dot_W2S = 0
            Q_dot_S2W = 0
            self.m[3] = self.m[4]
            self.T[4] = self.T[3]
            self.m[25] = self.m[4]
        
            self.h[3] = self.h[4]
            self.T[3] = temperature(P = self.P[3], H = self.h[3])

            # HPT 2
            self.h[5], h_HPT2, self.m[5], m_HPT2, x_HPT2, W_dot_HPT2 = turbine(m_HPT1, h_HPT1, self.P[3], self.P[5], eta_HPT, 0,0,0,0)
            self.T[5] = temperature(P = self.P[5], H = self.h[5])
            self.x[5] = quality(P = self.P[5], H = self.h[5])
            self.m[5] = m_HPT2
            self.h[5] = h_HPT2
            
            # IPT 1
            self.T[26] = self.T[23]+DT_HX_F
            self.P[26] = self.P[6]
            self.h[26] = enthalpy(T = self.T[26], P = self.P[26])
            self.h[6], h_IPT1, self.m[6], m_IPT1, x_IPT1, W_dot_IPT1 = turbine(self.m[5], self.h[5], self.P[5], self.P[6], eta_IPT1, Q_dot_HXF, self.h[26], self.m[25], self.h[25])
            self.m[26] = self.m[25]+self.m[6]
            self.T[6] = temperature(H = self.h[6], P = self.P[6])
            self.x[6] = quality(P = self.P[6], H = self.h[6])
        
            # IPT 2
            self.T[27] = self.T[22]+DT_HX_E
            self.P[27] = self.P[7]
            self.h[27] = enthalpy(T = self.T[27], P = self.P[27])
            self.h[7], h_IPT2, self.m[7], m_IPT2, x_IPT2, W_dot_IPT2 = turbine(m_IPT1, h_IPT1, self.P[6], self.P[7], eta_IPT2, Q_dot_HXE, self.h[27], self.m[26], self.h[26])
            self.m[27] = self.m[26]+self.m[7]
            self.T[7] = temperature(P = self.P[7], H = self.h[7])
            self.x[7] = quality(P = self.P[7], H = self.h[7])
        
            # IPT3
            self.P[21] = min(self.P[8],self.P[20],self.P[27])
            self.m[21] = self.m[1]
            self.m[8] = -(self.m[20] + self.m[27]) + self.m[21]
            self.h[8], W_dot_IPT3 = turbine_simple(m_IPT2, self.h[7], self.P[7], self.P[8], eta_HPT)
            self.h[21] = (self.m[20]*self.h[20] + (self.m[8]*self.h[8] + self.m[27]*self.h[27]))/self.m[21]
            self.h[22], W_dot_HP = pump(self.m[21], self.h[21], self.P[21], self.P[22], eta_IPT3)
            self.T[21] = temperature(P = self.P[21], H = self.h[21])
            self.T[8] = temperature(P = self.P[8], H = self.h[8])
            self.x[8] = quality(P = self.P[8], H = self.h[8])
            x_IPT3 = self.x[8]
        
            # LP Turbines     
            self.h[17] = enthalpy(P=self.P[17],T=self.T[17])
            Q_dot_HXA = self.m[16]*(self.h[17]-self.h[16])
            Q_dot_HXB = self.m[17]*(self.h[18]-self.h[17])
            Q_dot_HXC = self.m[18]*(self.h[19]-self.h[18])
            Q_dot_HXD = self.m[19]*(self.h[20]-self.h[19])

            # LPT 1
            self.h[9] = self.h[8]
            self.T[9] = self.T[8]
            self.m[9] = self.m[5] - self.m[6] - self.m[7] - self.m[8]
            self.T[28] = self.T[19] + DT_HX_D
            self.P[28] = self.P[10]
            self.h[28] = enthalpy(T = self.T[28], P = self.P[28])
            self.h[10], h_LPT1, self.m[10], m_LPT1, x_LPT1, W_dot_LPT1 = turbine(self.m[9], self.h[9], self.P[9], self.P[10], eta_LPT1, Q_dot_HXD, self.h[28], 0, 0)
            self.T[10] = temperature(P = self.P[10], H = self.h[10])
            self.x[10] = quality(P = self.P[10], H = self.h[10])
            self.m[28] = self.m[10]
        
            # LPT 2
            self.T[29] = self.T[18]+DT_HX_C
            self.P[29] = self.P[11]
            self.h[29] = enthalpy(T = self.T[29], P = self.P[29])
            self.h[11], h_LPT2, self.m[11], m_LPT2, x_LPT2, W_dot_LPT2 = turbine(m_LPT1, h_LPT1, self.P[10], self.P[11], eta_LPT2, Q_dot_HXC, self.h[29], self.m[28], self.h[28])
            self.T[11] = temperature(P = self.P[11], H = self.h[11])
            self.x[11] = quality(P = self.P[11], H = self.h[11])
            self.m[29] = self.m[28]+self.m[11]
        
            #LPT 3
            self.T[30] = self.T[17]+DT_HX_B
            self.P[30] = self.P[12]
            self.h[30] = enthalpy(T = self.T[30], P = self.P[30])
            self.h[12], h_LPT3, self.m[12], m_LPT3, x_LPT3, W_dot_LPT3 = turbine(m_LPT2, h_LPT2, self.P[11], self.P[12], eta_LPT3, Q_dot_HXB, self.h[30], self.m[29], self.h[29])
            self.T[12] = temperature(P = self.P[12], H = self.h[12])
            self.x[12] = quality(P = self.P[12], H = self.h[12])
            self.m[30] = self.m[29]+self.m[12]
        
            # LPT 4
            self.T[31] = self.T[16]+DT_HX_A
            self.P[31] = self.P[13]
            self.h[31] = enthalpy(T = self.T[31], P = self.P[31])
            self.h[13], h_LPT4, self.m[13], m_LPT4, x_LPT4, W_dot_LPT4 = turbine(m_LPT3, h_LPT3, self.P[12], self.P[13], eta_LPT4, Q_dot_HXA, self.h[31], self.m[30], self.h[30])
            self.T[13] = temperature(P = self.P[13], H = self.h[13])
            self.x[13] = quality(P = self.P[13], H = self.h[13])
            self.m[31] = self.m[30]+self.m[13]
        
            # LPT 5
            self.h[14],W_dot_LPT5 = turbine_simple(m_LPT4, h_LPT4, self.P[13], self.P[14], eta_LPT5)
            self.T[14] = temperature(P = self.P[14], H = self.h[14])
            self.x[14] = quality(P = self.P[14], H = self.h[14])
            self.m[14] = self.m[9]-self.m[10]-self.m[11]-self.m[12]-self.m[13]
            self.m[15] = self.m[14]+self.m[31]
        
            # condenser
            h_cond_in = (self.h[14]*self.m[14] + self.h[31]*self.m[31]) / (self.m[14]+self.m[31])
            Q_dot_cond = (self.m[14]+self.m[31])*(h_cond_in-self.h[15])

            # Efficiencies
            self.W_dot_HPT = W_dot_HPT1+W_dot_HPT2
            self.W_dot_IPT = W_dot_IPT1+W_dot_IPT2+W_dot_IPT3
            self.W_dot_LPT = W_dot_LPT1+W_dot_LPT2+W_dot_LPT3+W_dot_LPT4+W_dot_LPT5
            self.W_dot_gen = self.W_dot_HPT+self.W_dot_IPT+self.W_dot_LPT
            # Conservation of energy check
            eta_test = (self.W_dot_gen+Q_dot_cond+Q_dot_W2S)/(Q_dot_LFR+W_dot_LP+W_dot_HP+Q_dot_S2W)

            # --------------------------
            # //HX temperature plots
            # "HXF"
            UA_F_ds, UA_F_c, UA_F_dc, NTU_F, T_H_F, T_C_F = fwh(self.m[25], self.h[25], self.P[25], self.m[23], self.h[23], self.P[23], self.m[6],  self.h[6],  self.P[6],  DT_HX_F, self.N_hxrs)
            # "HXE"
            UA_E_ds, UA_E_c, UA_E_dc, NTU_E, T_H_E, T_C_E = fwh(self.m[26], self.h[26], self.P[26], self.m[22], self.h[22], self.P[22], self.m[7],  self.h[7],  self.P[7],  DT_HX_E, self.N_hxrs)
            # "HXD"
            UA_D_ds, UA_D_c, UA_D_dc, NTU_D, T_H_D, T_C_D = fwh(0,     0,     0,     self.m[19], self.h[19], self.P[19], self.m[10], self.h[10], self.P[10], DT_HX_D, self.N_hxrs)
            # "HXC"
            UA_C_ds, UA_C_c, UA_C_dc, NTU_C, T_H_C, T_C_C = fwh(self.m[28], self.h[28], self.P[28], self.m[18], self.h[18], self.P[18], self.m[11], self.h[11], self.P[11], DT_HX_C, self.N_hxrs)
            # "HXB"
            UA_B_ds, UA_B_c, UA_B_dc, NTU_B, T_H_B, T_C_B = fwh(self.m[29], self.h[29], self.P[29], self.m[17], self.h[17], self.P[17], self.m[12], self.h[12], self.P[12], DT_HX_B, self.N_hxrs)
            # "HXA"
            UA_A_ds, UA_A_c, UA_A_dc, NTU_A, T_H_A, T_C_A = fwh(self.m[30], self.h[30], self.P[30], self.m[16], self.h[16], self.P[16], self.m[13], self.h[13], self.P[13], DT_HX_A, self.N_hxrs)

            self.UA_table = {
                'A':[UA_A_ds,UA_A_c,UA_A_dc],
                'B':[UA_B_ds,UA_B_c,UA_B_dc],
                'C':[UA_C_ds,UA_C_c,UA_C_dc],
                'D':[UA_D_ds,UA_D_c,UA_D_dc],
                'E':[UA_E_ds,UA_E_c,UA_E_dc],
                'F':[UA_F_ds,UA_F_c,UA_F_dc],
            }
            for k in self.UA_table.keys():
                self.UA_table[k].append(sum(self.UA_table[k]))
            # --------------------------

            err_iter = 1 - eta_test

            # Compute all other temperature states
            for i in range(2,self.ns):
                if (self.h[i] != 0) and (self.P[i] != 0):
                    self.T[i] = temperature(H=self.h[i], P=self.P[i])

            #Efficiency ===============
            self.eta_cycle = self.W_dot_gen/Q_dot_LFR    # Cycle efficiency

            if self.logfile:
                print("    >> Time: {:5.1f}s | Error: {:.{}f} | Power: {:.1f} | eta: {:.4f} | ConsE: {:.4f}".format(time.time()-tstart, err_iter, abs(int(ceil(log10(tol))))+1, W_dot_gen, eta_cycle, eta_test))
                fout.write(self.__log_state(self.T,self.P,self.h,self.x,self.m))

            return err_iter
        # ----------------------

        root_scalar(cycle_eval, x0=self.m[15], x1=self.m[15]*1.01, maxiter=it_max, rtol=tol,)

        if self.logfile: fout.close()
    # --------------------
# ---------------

C = Cycle() #logfile='class-test-log.txt')
C.solve_cycle()
# r = fwh(C.m[30], C.h[30], C.P[30], C.m[16], C.h[16], C.P[16], C.m[13], C.h[13], C.P[13], 5.00, C.N_hxrs)
# fwh(C.m[25], C.h[25], C.P[25], C.m[23], C.h[23], C.P[23], C.m[6],  C.h[6],  C.P[6],  5.00, C.N_hxrs)
print(C)
asdf=1
