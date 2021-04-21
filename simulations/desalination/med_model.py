#MED Model for ReTI Research Lab under Professor Ben Lindley
#Created by Grace Stanke
import numpy as np
from iapws import SeaWater,IAPWS97
import pdb

class MED:
    def __init__(self):
        self.vapor_rate = [ 701.5, 700 ]                      
        self.brine_conc = .0335                              #Average concentration of salt in seawater
        self.brine_rate = 50.
        self.max_brine_conc= .067
        self.vapor_temp_n = [ 71.5, 70 ]                          #Known starting temp of the MED model - assuming we are starting at the second n-effect
        self.brine_temp_n = [ 23, 22 ]
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
        for i in range(1, 6):                    
            self.vapor_rate.append(self.vapor_flow_out(i))         
            self.brine_rate = self.brine_flow_out(i)                 #Updates brine_rate variable for every n-effect
            self.distill.append(sum(self.vapor_rate))                  #Adds the amount of distillate from each n-effect, creating a vector
        
        self.total_distill = sum(self.distill) #do this at the end of the loop
        
    def brine_flow_out(self,i):
        print(i)
        brine_out = (self.brine_rate*self.brine_conc)/(self.brine_rate-self.vapor_rate[i])
        return brine_out
        
    def vapor_flow_out(self, i):
        self.vapor_temp_n.append(self.vapor_temp_n[i-1] - self.tempchange)
        self.brine_temp_n.append(self.brine_temp_n[i-1] + self.tempchange)
        enth_bn = SeaWater(T=self.brine_temp_n[i]+273.15,S=self.brine_conc,P=self.pressure).h
        enth_vn = IAPWS97(T=self.vapor_temp_n[i]+273.15,P=self.pressure).h
        enth_bn_1 = SeaWater(T=self.brine_temp_n[i-1]+273.15,S=self.brine_conc,P=self.pressure).h
        enth_vn_1 = IAPWS97(T=self.vapor_temp_n[i-1]+273.15,P=self.pressure).h
        print(enth_bn, enth_vn, enth_bn_1, enth_vn_1)
        vapor_out = (1/(enth_vn - enth_bn))*(self.vapor_rate[i]*self.latentheat+(((self.distill[i]*self.max_brine_conc)\
            /(self.max_brine_conc-self.brine_conc))-self.distill[i])*(enth_bn_1- enth_bn))
        return vapor_out
        

    def feedflowrate(self, i):
        self.ffrate.append((self.max_brine_conc*self.vapor_rate[i])/(self.max_brine_conc-self.brine_conc))

MED()
