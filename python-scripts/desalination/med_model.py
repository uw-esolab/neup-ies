#MED Model for ReTI Research Lab under Professor Ben Lindley
#Created by Grace Stanke
import numpy as np
from iapws import SeaWater,IAPWS97

"""Things to fix!
Currently get a divide by zero in brine_flow_out

You appear to be cycling through temperature in vapor_flow_out. Set temperature
as a list i.e. self.temp = [] in __init__
There are two ways of handling this in vapor_flow_out. You need to pick one.
EITHER use the class variable in the function and 
pass it i (as now done in brine_flow_out) so you have things like self.temp[i] 
and self.temp[i+1] 
OR you can have the function work with LOCAL variables temp_n and temp_n_1 -
if you do this, then you need to return the next ste of temperatures to the line
that calls it, and assign them to the right variable there.

It's probably less work to use the first method


"""

class MED:
    def __init__(self):
        self.vapor_rate = [ 100., 90 ]                      
        self.brine_conc = .035                              #Average concentration of salt in seawater
        self.brine_rate = 100.
        self.max_brine_conc=100.
        self.temp_n_1 = 37                                  #Known starting temp of the MED model - 37 C is assuming we are starting at the second n-effect
        self.latentheat = 100.
        self.tempchange = 3.                                #Known temperature change from Sharan paper
        self.distill = [ 0, 1 ]                             #starting value
        self.pressure = 0.101325 #TODO: check this is in MPa
        self.salinity = 0.001 #TODO: find actual value
        # self.vapor_rate = vapor flow rate
        # self.brine_conc = brine concentration
        # self.max_brine_conc = maximum brine concentration
        # self.brine_rate = brine flow rate
        # self.temp_n_1 = starting temperature, assume to be 40 C
        # self.latentheat = latent heat 
        # self.tempchange = temperature change between each n effect
        for i in range(5):                                           #this is more standard
            self.brine_rate += self.brine_flow_out(i)                 #Updates brine_rate variable for every n-effect
            
            #BL: you had a function called the same thing as a variable so I've put it as one line
            self.distill.append(sum(self.vapor_rate))                  #Adds the amount of distillate from each n-effect, creating a vector
            self.vapor_rate.append(self.vapor_flow_out())       #Adds the vapor rate from each effect to the vector
        
        self.total_distill = sum(self.distill) #do this at the end fo the loop
        
    def brine_flow_out(self,i):
        brine_out = (self.brine_rate*self.brine_conc)/(self.brine_rate-self.vapor_rate[i])
        return brine_out
        
    def vapor_flow_out(self):
        temp_n = temp_n_1 - tempchange
        #CHECK: does pressure change?
        enth_bn = SeaWater(T=temp_n,S=self.salinity,P=self.pressure) 
        enth_vn = IAPWS97(T=temp_n,P=self.pressure)
        enth_bn_1 = SeaWater(T=temp_n_1,S=self.salinity,P=self.pressure)
        enth_vn_1 = IAPWS97(T=temp_n_1,P=self.pressure)
        sumvpfr = np.sum(vapor_rate)
        vapor_flow_out = (1/(enth_vn-enth_bn))*(vapor_rate[i]*latentheat+(((distill*max_brine_conc)/(max_brine_conc-feed_conc))-sumvpfr)*(enth_bn_1-enth_bn))
        self.temp_n_1 += temp_n

    def feedflowrate(self):
            ffrate[i] = (max_brine_conce*vapor_rate[i])/(max_brine_conc-brine_conc)


