"""
MED Model for ReTI Research Lab under Professor Ben Lindley
Created by Grace Stanke and Ben Lindley
Based on Sharan 2019
As with Sharan, gives optimal water output with 4 effects
"""

import numpy as np
from iapws import SeaWater,IAPWS97

def water_value_in_arizona():
    """https://efc.web.unc.edu/2014/09/23/conservation-water-rates-arizona-utilities-using-rates-discourage-discretionary-water-use/"""
    price_per_1000_gallons=3.5
    gallons_per_m3 = 264.172
    price_per_m3 = price_per_1000_gallons/1000*gallons_per_m3
    return price_per_m3

class MED:
    
    """values we ultimately expect to pass to the MED model:
        k=4
        water rate - solved through EES mass balance
        water_temp - solved through EES mass and temperature balance
    """
    
    def __init__(self,k):
        self.vapor_rate = []                          #Flow of the vapor rate for each n-effect
        self.feed_conc = .0335                       #Average concentration of salt in seawater
        self.brine_rate = []                         #Brine flow rate 
        self.max_brine_conc= .067                     #Maximum brine concentration
        self.k = k                                   #Number of desired effects + 1 for starting values
        self.vbtemp = np.zeros(self.k)            #Vapor and brine temperature vector creation
        self.vbtemp[self.k-1] = 40                  #Fixing the final temperature of the brine and vapor
        self.enth_vapor = np.zeros(self.k)            #Vapor enthalpy vector creation
        self.enth_brine = np.zeros(self.k)            #Brine enthalpy vector creation
        self.enth_brine_nea = np.zeros(self.k)       #allowance for flashing
        self.water_temp = 67.6                     #Known starting/inlet temp of the DEMINERALIZED WATER [T2] (sCO2 outlet - delta_T_PCHE)
        #self.feed_temp  = 20                       #Known starting/inlet temp of the BRINE [T1]
        self.latentheat = np.zeros(self.k)            #Latent heat of vapor vector creation
        self.tempchange = 3.                          #Known temperature change, from Sharan paper (delta_T_NEA)
        self.water_rate = 308                      #Is this the feed flow rate of the DEMINERALIZED WATER? m_dot_w2
        self.pressure = np.zeros(self.k)                          #Pressure in MPa 
        self.brine_out = []
        self.nea = 3. 
        """NEA: upper value from MSF NEA is 1. The correlation Sharan used (JPME paper) is very dififcult to interpret. 
        te difference is fairly small, perhaps 4% higher than Sharan with 4 effects. Setting NEA to 3 gives
        better agreement and could be due to looking at TD-TF in this paper, which is about 3C
        this gives almost perfect agreement with Sharan for 4 effects"""
        
        brine_conc_in = [self.feed_conc for i in range(k)]
        
        max_dif=1
        """we had to assume brine concentrations to get the enthalpies. We iterate the brine concentraitons (one loop does it)
        until we get convergence"""
        while max_dif>0.01:
            brine_conc_out = self.iterate_MED(k,brine_conc_in) #calculate Brine concentartion
            max_dif=np.max(np.abs(np.divide(np.array(brine_conc_out),np.array(brine_conc_in))-1)) #difference to previous iteration
            brine_conc_in=brine_conc_out
            
         
    
    def iterate_MED(self,k,brine_conc_in):

        #step backwards through the stages to get the vapor temperatures
        for i in range(self.k-2, -1, -1):
            self.vbtemp[i] = (self.vbtemp[i+1] + self.tempchange)
        
        self.feed_temp = self.vbtemp[0] # A closer read of Sharan suggests this is the temperature of the first effect
        
        
        for i in range(0, self.k):
            
            self.enth_vapor[i] = IAPWS97(T=self.vbtemp[i]+273.15,x=1).h
            self.pressure[i] = IAPWS97(T=self.vbtemp[i]+273.15,x=1).P
            #self.enth_brine[i] = SeaWater(T=self.vbtemp[i]+273.15,S=self.feed_conc,P=self.pressure[i]).h
            self.enth_brine[i] = self.seaWaterSatH(self.vbtemp[i]+273.15,brine_conc_in[i],self.pressure[i])
            self.enth_brine_nea[i] = self.seaWaterSatH(self.vbtemp[i]+273.15-self.nea,brine_conc_in[i],self.pressure[i])
            self.latentheat[i] = self.enth_vapor[i]-IAPWS97(T=self.vbtemp[i]+273.15,x=0).h
            
        self.enth_feed =  self.seaWaterSatH(self.feed_temp+273.15,self.feed_conc,self.pressure[i])
        

        
        ##Finding the vapor_rate for each n-effect. See ReTI/NE-2/Technical/Desalination Equations from Sharan 
        #for definitions of A and C matrices
        C = np.zeros(self.k)
        A = np.zeros((self.k,self.k))
        A[0,0] = (-self.max_brine_conc/(self.max_brine_conc-self.feed_conc))*\
                        (self.enth_feed - self.enth_brine[0])+(self.enth_vapor[0]-self.enth_brine[0])
        
        self.cp_water = IAPWS97(T=self.vbtemp[0]+273.15-0.1,P=self.pressure[0]).cp
        C[0] = self.water_rate*(self.water_temp-(self.feed_temp+self.tempchange))*self.cp_water
        
        for j in range(1,self.k):           #Creating the first row of the matrix
            A[0,j] = (-self.max_brine_conc/(self.max_brine_conc-self.feed_conc))*\
                        (self.enth_feed - self.enth_brine[0])
                        
        for q in range(1,self.k):           #Filling in the rest of the matrix
            #Below is the diagonal of the matrix
            A[q,q]= -(self.enth_vapor[q] - self.enth_brine[q]) + \
                (self.max_brine_conc/(self.max_brine_conc-self.feed_conc))*\
                    (self.enth_brine_nea[q-1] - self.enth_brine[q])
            #Below is on column to the left of the diagonal of the matrix
            A[q,q-1]= self.latentheat[q-1]+((self.max_brine_conc/(self.max_brine_conc-self.feed_conc))-1)*\
                (self.enth_brine_nea[q-1] - self.enth_brine[q])
            for j in range(0,q-1):
                #Below is everything to the left, under the q-1 diagonal. 
                A[q,j]= (self.max_brine_conc/(self.max_brine_conc-self.feed_conc)-1)*\
                        (self.enth_brine_nea[q-1] - self.enth_brine[q])
            for j in range(q+1,self.k):
                #Below is everything to the right of the q diagonal
                A[q,j]= (self.max_brine_conc/(self.max_brine_conc-self.feed_conc))*\
                        (self.enth_brine_nea[q-1] - self.enth_brine[q])
                        
            #For eqn 2-n, every term is related to a vapor rate
            C[q] = 0
            
        invA = np.linalg.inv(A) 
        self.vapor_rate = np.dot(invA,C) #Solve for the vapor rate
        self.distill = sum(self.vapor_rate)
        
        #Calculate Brine feed flow with Eq 10
        self.feed_rate = self.distill*self.max_brine_conc/(self.max_brine_conc-self.feed_conc)

        brine_conc = []
        for i in range(self.k):                          
            self.brine_rate.append(self.brine_flow_out(i))    #Updates brine_rate variable for every n-effect
            brine_conc.append(self.brine_conctn(i))
        
        self.m3_per_day = self.distill*3600*24/1000.0 #from Sharan - convert the units from m3/s
        return brine_conc
        
    
        
    def brine_flow_out(self,i):
        brine_out = self.feed_rate-sum(self.vapor_rate[:i+1]) #Eq 6
        return brine_out
    
    def brine_conctn(self,i):
        brine_concentration = (self.feed_rate*self.feed_conc)/(self.feed_rate-sum(self.vapor_rate[:i+1])) #Eq 8
        return brine_concentration
       
    def seaWaterSatH(self,T,S,P):
        """the enthalpy at saturation pressure was pinging between gas and liquid. This clears it up"""
        test_enth = SeaWater(T=T,P=P,S=S).h
        if test_enth < 1000:
            return test_enth
        else:
            return self.seaWaterSatH(T-0.1,S,P)

for k in range(1,10):
    x=MED(k)
    print(x.m3_per_day)



    
