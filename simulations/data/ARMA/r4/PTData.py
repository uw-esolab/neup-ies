# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:55:45 2023

@author: aidan
"""
from PTPlant import ptrun
import matplotlib.pyplot as plt

synfiles = []
realfiles=[]
variable = 'PT'

annualenergy_list_s = []
capfactor_list_s = []
ppa_list_s = []
lcoe_list_s = []
npv_list_s = []
cost_list_s = []    
annualenergy_list_r = []
capfactor_list_r = []
ppa_list_r = []
lcoe_list_r = []
npv_list_r = []
cost_list_r = []    

for i in range(1,100,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/TestSyn' + str(i) + '.csv'
    fullfile = bytes(fullfile,'utf-8')
    synfiles.append(fullfile)
    
for i in range(1998,2020,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/Weather_Data/102574_35.93_-115.26_' + str(i) + '.csv' 
    fullfile = bytes(fullfile,'utf-8')
    realfiles.append(fullfile)
   
for file in synfiles:    
    annualenergy, capfactor, ppa, lcoe, npv, cost = ptrun(file)
    
    annualenergy_list_s.append(annualenergy)
    capfactor_list_s.append(capfactor)
    ppa_list_s.append(ppa)
    lcoe_list_s.append(lcoe)
    npv_list_s.append(npv)
    cost_list_s.append(cost)
    
for file in realfiles:    
    annualenergy, capfactor, ppa, lcoe, npv, cost = ptrun(file)
    
    annualenergy_list_r.append(annualenergy)
    capfactor_list_r.append(capfactor)
    ppa_list_r.append(ppa)
    lcoe_list_r.append(lcoe)
    npv_list_r.append(npv)
    cost_list_r.append(cost)
    
print(annualenergy_list_s)
print(capfactor_list_s)
print(ppa_list_s)
print(lcoe_list_s)
print(npv_list_s)
print(cost_list_s)    
print(annualenergy_list_r)
print(capfactor_list_r)
print(ppa_list_r)
print(lcoe_list_r)
print(npv_list_r)
print(cost_list_r)

plt.hist(annualenergy_list_s, bins=10, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
plt.hist(annualenergy_list_r, bins=10, facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
plt.title("Histogram for each year of Annual Energy Production")
plt.savefig('Figures/Annual Energy Production'+variable+'.png')
plt.ylabel('Frequency')
plt.xlabel('Annual Energy Production / kWh')
plt.legend()
plt.show()

plt.hist(capfactor_list_s, bins=10, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
plt.hist(capfactor_list_r, bins=10, facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
plt.title("Histogram for each year of Capacity Factor")
plt.ylabel('Frequency')
plt.xlabel('Capacity Factor / %')
plt.savefig('Figures/CapacityFactor'+variable+'.png')
plt.legend()
plt.show()

plt.hist(ppa_list_s, bins=10, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
plt.hist(ppa_list_r, bins=10, facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
plt.title("Histogram for each year of PPA")
plt.ylabel('Frequency')
plt.xlabel('PPA / c/kWh')
plt.savefig('Figures/PPA'+variable+'.png')
plt.legend()
plt.show()

plt.hist(lcoe_list_s, bins=10, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
plt.hist(lcoe_list_r, bins=10, facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
plt.title("Histogram for each year of LCOE")
plt.ylabel('Frequency')
plt.xlabel('LCOE /  c/kWh')
plt.savefig('Figures/LCOE'+variable+'.png')
plt.legend()
plt.show()

plt.hist(npv_list_s, bins=10, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
plt.hist(npv_list_r, bins=10, facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
plt.title("Histogram for each year of NPV")
plt.ylabel('Frequency')
plt.xlabel('NPV / $ ')
plt.savefig('Figures/NPV'+variable+'.png')
plt.legend()
plt.show()

plt.hist(cost_list_s, bins=10, facecolor = 'none', edgecolor='blue', label = 'Synthetic data', density = True)
plt.hist(cost_list_r, bins=10, facecolor = 'none', edgecolor='red', label = 'Real data', density = True)
plt.title("Histogram for each year of Net Capital Cost")
plt.ylabel('Frequency')
plt.xlabel('Net Capital Cost / $ ')
plt.savefig('Figures/NCC'+variable+'.png')
plt.legend()
plt.show()