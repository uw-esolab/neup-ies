# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 08:59:45 2026

@author: b9801
"""

import zld_model_v3 as zld
import numpy as np
import matplotlib.pyplot as plt

gather_data1 = False
gather_data2 = False
gather_data3 = False
Li_prices = 50*np.arange(20)
ED_elec_reqs = 0.5*np.arange(16)
RO_elec_reqs = 2.5+0.25*np.arange(9)

if gather_data1:

    Li = np.zeros((Li_prices.shape[0],ED_elec_reqs.shape[0]))
    MED = np.zeros((Li_prices.shape[0],ED_elec_reqs.shape[0]))
    RO = np.zeros((Li_prices.shape[0],ED_elec_reqs.shape[0]))
    
    for j,Li_price in enumerate(Li_prices):
        for k,ED_elec_req in enumerate(ED_elec_reqs):
            Li[j,k],MED[j,k],RO[j,k]=zld.main(Li_price=Li_price,ED_elec_req=ED_elec_req,RO_elec_req=45) #force RO out of system so stable answer
        print(j)
        
    np.savez("tmp_results.npz", Li=Li,MED=MED,RO=RO)
    
else:
    np.load("tmp_results.npz")
    
fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)

pcm1 = axs[0].pcolormesh(ED_elec_reqs, Li_prices,  Li, shading='auto')
fig.colorbar(pcm1, ax=axs[0])
axs[0].set_title("Li Recovered (kg/hr)")
axs[0].set_ylabel("Li Yield x Price ($/kg)")
axs[0].set_xlabel("ED Electricity Requirement (kWh/m³)")

pcm2 = axs[1].pcolormesh(ED_elec_reqs, Li_prices, MED, shading='auto')
fig.colorbar(pcm2, ax=axs[1])
axs[1].set_title("MED Capacity (m³/h)")
axs[1].set_xlabel("ED Electricity Requirement (kWh/m³)")

plt.tight_layout()
plt.show()
fig.savefig("physor_updated.png",dpi=1000)
    
if gather_data2:
        
    Li2 = np.zeros((RO_elec_reqs.shape[0],ED_elec_reqs.shape[0]))
    MED2 = np.zeros((RO_elec_reqs.shape[0],ED_elec_reqs.shape[0]))
    RO2 = np.zeros((RO_elec_reqs.shape[0],ED_elec_reqs.shape[0]))
    
    for j,RO_elec_req in enumerate(RO_elec_reqs):
        for k,ED_elec_req in enumerate(ED_elec_reqs):
            Li2[j,k],MED2[j,k],RO2[j,k]=zld.main(Li_price=700,ED_elec_req=ED_elec_req,RO_elec_req=RO_elec_req) #force RO out of system so stable answer
        print(j)
    
    np.savez("tmp_results2.npz",Li2=Li2,MED2=MED2,RO2=RO2)
    
else:
    np.load("tmp_results2.npz")
    
fig, axs = plt.subplots(1, 3, figsize=(10, 4), sharex=True, sharey=True)

pcm1 = axs[0].pcolormesh(ED_elec_reqs, RO_elec_reqs,  Li2, shading='auto')
fig.colorbar(pcm1, ax=axs[0])
axs[0].set_title("Li Recovered (kg/hr)")
axs[0].set_ylabel("RO Electricity Requirement (kWh/m³)")
axs[0].set_xlabel("ED Electricity Requirement (kWh/m³)")

pcm2 = axs[1].pcolormesh(ED_elec_reqs, RO_elec_reqs, MED2, shading='auto')
fig.colorbar(pcm2, ax=axs[1])
axs[1].set_title("MED Capacity (m³/h)")
axs[1].set_xlabel("ED Electricity Requirement (kWh/m³)")

pcm3 = axs[2].pcolormesh(ED_elec_reqs, RO_elec_reqs, RO2, shading='auto')
fig.colorbar(pcm3, ax=axs[2], label="RO Capacity (m³/h)")
axs[2].set_title("RO Capacity (m³/h)")
axs[2].set_xlabel("ED Electricity Requirement (kWh/m³)")

plt.tight_layout()
plt.show()
fig.savefig("physor_new.png",dpi=1000)

    
if gather_data3:
    
    Li3 = np.zeros((Li_prices.shape[0],ED_elec_reqs.shape[0]))
    MED3 = np.zeros((Li_prices.shape[0],ED_elec_reqs.shape[0]))
    RO3 = np.zeros((Li_prices.shape[0],ED_elec_reqs.shape[0])) 
    
    for j,Li_price in enumerate(Li_prices):
        for k,ED_elec_req in enumerate(ED_elec_reqs):
            Li3[j,k],MED3[j,k],RO3[j,k]=zld.main(Li_price=Li_price,ED_elec_req=ED_elec_req,RO_elec_req=2.0)#,MED_elec_req=10.0) #force RO out of system so stable answer
        print(j)
        
    np.savez("tmp_results3.npz", Li3=Li3,MED3=MED3,RO3=RO3)
    
else:
    np.load("tmp_results3.npz")



fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)

pcm1 = axs[0].pcolormesh(ED_elec_reqs, Li_prices,  Li3, shading='auto')
fig.colorbar(pcm1, ax=axs[0])
axs[0].set_title("Li Recovered (kg/hr)")
axs[0].set_ylabel("Li Yield x Price ($/kg)")
axs[0].set_xlabel("ED Electricity Requirement (kWh/m³)")

pcm2 = axs[1].pcolormesh(ED_elec_reqs, Li_prices, RO3, shading='auto')
fig.colorbar(pcm2, ax=axs[1])
axs[1].set_title("RO Capacity (m³/h)")
axs[1].set_xlabel("ED Electricity Requirement (kWh/m³)")

plt.tight_layout()
plt.show()
fig.savefig("physor_new2.png",dpi=1000)

