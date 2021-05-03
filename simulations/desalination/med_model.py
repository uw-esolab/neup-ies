#MED Model for ReTI Research Lab under Professor Ben Lindley
#Created by Grace Stanke
import numpy as np
from iapws import SeaWater,IAPWS97
from numpy.linalg import inv

class MED:
    def __init__(self):
        self.vapor_rate = [ 701.5 ]                   #Flow of the vapor rate for each n-effect
        self.brine_conc = .0335                       #Average concentration of salt in seawater
        self.brine_rate = 50.                         #Brine flow rate, constant through each n effect
        self.max_brine_conc= .067                     #Maximum brine concentration
        self.k = 6                                    #Number of desired effects + 1 for starting values
        self.vapor_temp = np.zeros(self.k)            #Vapor temperature vector creation
        self.brine_temp = np.zeros(self.k)            #Brine temparature vector creation
        self.enth_vapor = np.zeros(self.k)            #Vapor enthalpy vector creation
        self.enth_brine = np.zeros(self.k)            #Brine enthalpy vector creation
        self.vapor_temp[0] = 71.5                     #Known starting/inlet temp of the vapor in C
        self.brine_temp[0] = 23                       #Known starting/inlet temp of the brine in C
        self.latentheat = np.zeros(self.k)            #Latent heat of vapor vector creation
        self.latentheat[0] = 2446.4                   #Starting latent heat value at Brine temp (0)
        self.tempchange = 3.                          #Known temperature change, from Sharan paper
        self.ffrate = [ ]                             #Starting feed flow rate - NOT NEEDED FOR THIS PROJECT
        self.pressure = 7.75                          #Pressure in MPa
        
        for i in range(1, self.k):
            self.vapor_temp[i] = (self.vapor_temp[i-1] - self.tempchange)
            self.brine_temp[i] = (self.brine_temp[i-1] + self.tempchange)
            self.enth_brine[i] = SeaWater(T=self.brine_temp[i]+273.15,S=self.brine_conc,P=self.pressure).h
            self.enth_vapor[i] = IAPWS97(T=self.vapor_temp[i]+273.15,P=self.pressure).h
            self.latentheat[i] = -2.36985*self.brine_temp[i] + 2500.9
        
        
        ##Finding the vapor_rate for each n-effect
        C = np.zeros(self.k)
        A = np.zeros((self.k,self.k))
        A[0,0] = (-self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (self.enth_brine[0] - self.enth_brine[1])
        C[0] = self.vapor_rate[0]*(self.latentheat[0] - self.enth_vapor[0]+\
                        self.enth_brine[0]+((-self.max_brine_conc/(self.max_brine_conc-self.brine_conc))-1)*\
                        (self.enth_brine[0]-self.enth_brine[1]))
        
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
            #Below is the construction of the C vector. 
            C[q] = self.vapor_rate[0]*(self.max_brine_conc/(self.max_brine_conc-self.brine_conc)-1)*\
                        (self.enth_brine[q-1] - self.enth_brine[q])
        invA = inv(A)
        self.vapor_rate = np.dot(invA,C)
        
        for i in range(1, self.k):                          
            self.brine_rate = self.brine_flow_out(i)    #Updates brine_rate variable for every n-effect
            self.distill = sum(self.vapor_rate)
        
    def brine_flow_out(self,i):
        brine_out = (self.brine_rate*self.brine_conc)/(self.brine_rate-self.vapor_rate[i])
        return brine_out

    def feedflowrate(self, i):          #Unused function for now, but leaving it to simply exist for this moment. 
        self.ffrate.append((self.max_brine_conc*self.vapor_rate[i])/(self.max_brine_conc-self.brine_conc))

x = MED()
print(x.distill)
