#MED Model for ReTI Research Lab under Professor Ben Lindley
#Created by Grace Stanke
import numpy as np
from iapws import SeaWater,IAPWS97
import pdb

class MED:
    def __init__(self):
        self.vapor_rate = [ 701.5 ]                      
        self.brine_conc = .0335                              #Average concentration of salt in seawater
        self.brine_rate = 50.
        self.max_brine_conc= .067
        self.k = 6                                          #Number of desired effects + 1
        self.vapor_temp = np.zeros(k)
        self.brine_temp = np.zeros(k)
        self.enth_vapor = np.zeros(k)
        self.enth_brine = np.zeros(k)
        self.vapor_temp[0] = 71.5                          #Known starting temp of the MED model - assuming we are starting at the second n-effect
        self.brine_temp[0] = 23
        self.latentheat = 10.
        self.tempchange = 3.                                #Known temperature change from Sharan paper
        self.distill = [ 0, 1 ]
        self.ffrate = [ ]                             #starting value
        self.pressure = 7.75 #TODO: check this is in MPa
        # self.vapor_rate = vapor flow rate
        # self.brine_conc = brine concentration
        # self.max_brine_conc = maximum brine concentration
        # self.brine_rate = brine flow rate
        # self.temp_n_1 = starting temperature, assume to be 40 C
        # self.latentheat = latent heat 
        # self.tempchange = temperature change between each n effect
        for i in range(1, self.k):                    
            self.vapor_rate.append(self.vapor_flow_out(i))         
            self.brine_rate = self.brine_flow_out(i)                 #Updates brine_rate variable for every n-effect
            self.distill.append(sum(self.vapor_rate))                  #Adds the amount of distillate from each n-effect, creating a vector
        
        self.total_distill = sum(self.distill) #do this at the end of the loop
        
    def brine_flow_out(self,i):
        print(i)
        brine_out = (self.brine_rate*self.brine_conc)/(self.brine_rate-self.vapor_rate[i])
        return brine_out
        
    def vapor_flow_out(self, i):
        self.vapor_temp[i] = (self.vapor_temp[i-1] - self.tempchange)
        self.brine_temp[i] = (self.brine_temp[i-1] + self.tempchange)
        self.enth_brine[i] = SeaWater(T=self.brine_temp[i]+273.15,S=self.brine_conc,P=self.pressure).h
        self.enth_vapor[i] = IAPWS97(T=self.vapor_temp[i]+273.15,P=self.pressure).h
        #vapor_out = (1/(enth_vn - enth_bn))*(self.vapor_rate[i]*self.latentheat+(((self.distill[i]*self.max_brine_conc)\
        #   /(self.max_brine_conc-self.brine_conc))-self.distill[i])*(enth_bn_1- enth_bn))
        C = np.zeros(self.K)
        A = np.zeros(self.K,self.K)
        A[0,0] = (-self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (enth_brine[0] - enth_brine[1])
        C[0] = <expression for C_1>
        #complete first row
        for j in range(1,self.K):
        A[0,j] = (-self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (enth_brine[0] - enth_brine[1])
        #now do all the other rows
        for q in range(1,self.K):
            A[q,q]= -(enth_vapor[q] - enth_brine[q]) + \
                (self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                    (enth_brine[q-1] - enth_brine[q])
            A[q,q-1]= <expression for A_n,(n-1)>
            for j in range(0,q-1):
                A[q,j]= (self.max_brine_conc/(self.max_brine_conc-self.brine_conc)-1)*\
                        (enth_brine[q-1] - enth_brine[q])
                for j in range(q+1,self.K):
                    A[q,j]= (self.max_brine_conc/(self.max_brine_conc-self.brine_conc))*\
                        (enth_brine[q-1] - enth_brine[q])
                    C[q] = self.vapor_rate[0]*(self.max_brine_conc/(self.max_brine_conc-self.brine_conc)-1)\
                        (enth_brine[q-1] - enth_brine[q])
                    self.vapor_rate = np.dot(np.inv(A),C)
        return vapor_out
        

    def feedflowrate(self, i):
        self.ffrate.append((self.max_brine_conc*self.vapor_rate[i])/(self.max_brine_conc-self.brine_conc))

MED()
