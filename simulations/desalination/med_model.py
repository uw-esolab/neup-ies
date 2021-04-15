#MED Model for ReTI Research Lab under Professor Ben Lindley
#Created by Grace Stanke
import numpy as np
from iapws import SeaWater,IAPWS97

class MED:
    def __init__(self):
        self.vapor_rate = [ 701.5 ]                      
        self.brine_conc = .0335                              #Average concentration of salt in seawater
        self.brine_rate = 50.
        self.max_brine_conc= .067
        self.vapor_temp_n = [ 71.5 ]                          #Known starting temp of the MED model - assuming we are starting at the second n-effect
        self.brine_temp_n = [ 23 ]
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
        for i in range(5):                                           #this is more standard
            self.brine_rate = self.brine_flow_out(i)                 #Updates brine_rate variable for every n-effect
            self.vapor_rate.append(self.vapor_flow_out(i))       #Adds the vapor rate from each effect to the vector
            self.distill.append(sum(self.vapor_rate))                  #Adds the amount of distillate from each n-effect, creating a vector
        
        self.total_distill = sum(self.distill) #do this at the end of the loop
        
    def brine_flow_out(self,i):
        brine_out = (self.brine_rate*self.brine_conc)/(self.brine_rate-self.vapor_rate[i])
        return brine_out
        
    def vapor_flow_out(self, i):
        self.vapor_temp_n[i] = self.vapor_temp_n[i-1] - self.tempchange
        self.brine_temp_n[i] = self.brine_temp_n[i-1] + self.tempchange
        enth_bn = SeaWater(T=self.brine_temp_n[i]+273.15,S=self.brine_conc,P=self.pressure)
        enth_vn = IAPWS97(T=self.vapor_temp_n[i]+273.15,P=self.pressure)
        enth_bn_1 = SeaWater(T=self.brine_temp_n[i-1]+273.15,S=self.brine_conc,P=self.pressure)
        enth_vn_1 = IAPWS97(T=self.vapor_temp_n[i-1]+273.15,P=self.pressure)
        return enth_bn
        vapor_out = (1/(np.asarray(enth_vn.h) - np.asarray(enth_bn.h)))*(self.vapor_rate[i]*self.latentheat+(((self.distill*self.max_brine_conc)/(self.max_brine_conc-self.brine_conc))-self.distill)*(np.asarray(enth_bn_1.h)- np.asarray(enth_bn.h)))
        return vapor_out
        

    def feedflowrate(self, i):
        self.ffrate.append((self.max_brine_conc*self.vapor_rate[i])/(self.max_brine_conc-self.brine_conc))

MED()
