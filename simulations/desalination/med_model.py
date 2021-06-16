#MED Model for ReTI Research Lab under Professor Ben Lindley
#Created by Grace Stanke
import numpy as np
from iapws import SeaWater,IAPWS97

class MED:
    
    def __init__(self):
        self.vapor_rate = []                          #Flow of the vapor rate for each n-effect
        self.brine_conc = .0335                       #Average concentration of salt in seawater
        self.brine_rate = []                         #Brine flow rate 
        self.max_brine_conc= .067                     #Maximum brine concentration
        self.k = 6                                    #Number of desired effects + 1 for starting values
        self.vbtemp = np.zeros(self.k)            #Vapor and brine temperature vector creation
        self.vbtemp[self.k-1] = 40                  #Fixing the final temperature of the brine and vapor
        self.enth_vapor = np.zeros(self.k)            #Vapor enthalpy vector creation
        self.enth_brine = np.zeros(self.k)            #Brine enthalpy vector creation
        self.water_temp = 68.5                     #Known starting/inlet temp of the DEMINERALIZED WATER [T2] (sCO2 outlet - delta_T_PCHE)
        self.feed_temp  = 20                       #Known starting/inlet temp of the BRINE [T1]
        self.latentheat = np.zeros(self.k)            #Latent heat of vapor vector creation
        self.tempchange = 3.                          #Known temperature change, from Sharan paper (delta_T_NEA)
        self.water_rate = 701.5                        #Is this the feed flow rate of the DEMINERALIZED WATER? m_dot_w2
        self.pressure = np.zeros(self.k)                          #Pressure in MPa 
        
        
        for i in range(self.k-2, -1, -1):
            self.vbtemp[i] = (self.vbtemp[i+1] + self.tempchange)
            
        
        for i in range(0, self.k):

            self.enth_vapor[i] = IAPWS97(T=self.vbtemp[i]+273.15,x=1).h
            self.pressure[i] = IAPWS97(T=self.vbtemp[i]+273.15,x=1).P
            self.enth_brine[i] = SeaWater(T=self.vbtemp[i]+273.15,S=self.brine_conc,P=self.pressure[i]).h
            self.latentheat[i] = -2.36985*self.vbtemp[i] + 2500.9       #Equation found from a linear relation in EES
        self.enth_feed =  SeaWater(T=self.feed_temp+273.15,S=self.brine_conc,P=self.pressure[0]).h
        
        ##Finding the vapor_rate for each n-effect
        C = np.zeros(self.k)
        A = np.zeros((self.k,self.k))
        A[0,0] = (-self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (self.enth_feed - self.enth_brine[0])+(self.enth_vapor[0]-self.enth_brine[0])
        
        self.cp_water = IAPWS97(T=self.vbtemp[0]+273.15,P=self.pressure[0]).cp
        C[0] = self.water_rate*(self.water_temp-(self.feed_temp+self.tempchange))*self.cp_water
        
        for j in range(1,self.k):           #Creating the first row of the matrix
            A[0,j] = (-self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (self.enth_brine[0] - self.enth_brine[1])
                        
        for q in range(1,self.k):           #Filling in the rest of the matrix
            #Below is the diagonal of the matrix
            A[q,q]= -(self.enth_vapor[q] - self.enth_brine[q]) + \
                (self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                    (self.enth_brine[q-1] - self.enth_brine[q])
            #Below is on column to the left of the diagonal of the matrix
            A[q,q-1]= self.latentheat[q-1]+((self.max_brine_conc/(self.max_brine_conc-self.brine_conc))-1)*\
                (self.enth_brine[q-1] - self.enth_brine[q])
            for j in range(0,q-1):
                #Below is everything to the left, under the q-1 diagonal. 
                A[q,j]= (self.max_brine_conc/(self.max_brine_conc-self.brine_conc)-1)*\
                        (self.enth_brine[q-1] - self.enth_brine[q])
            for j in range(q+1,self.k):
                #Below is everything to the right of the q diagonal
                A[q,j]= (self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (self.enth_brine[q-1] - self.enth_brine[q])
                        
            #For eqn 2-n, every term is related to a vapor rate
            C[q] = 0
            
        invA = np.linalg.inv(A) 
        self.vapor_rate = np.dot(invA,C)
        self.distill = sum(self.vapor_rate)
        
        #Calculate Brine feed flow with Eq 10
        self.feed_rate = self.distill*self.max_brine_conc/(self.max_brine_conc-self.brine_conc)
        for i in range(self.k):                          
            self.brine_rate.append(self.brine_flow_out(i))    #Updates brine_rate variable for every n-effect
        
        
    def brine_flow_out(self,i):
        brine_out = self.feed_rate-sum(self.vapor_rate[:i+1]) #Eq 6
        return brine_out

    def feedflowrate(self, i):          #Unused function for now, but leaving it to simply exist for this moment. 
        self.ffrate.append((self.max_brine_conc*self.vapor_rate[i])/(self.max_brine_conc-self.brine_conc))

x = MED()
print(x.distill)
