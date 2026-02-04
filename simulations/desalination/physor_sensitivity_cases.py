# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 08:59:45 2026

@author: b9801
"""

import zld_model_v3 as zld
import numpy as np
import matplotlib.pyplot as plt

gather_data = False
Li_prices = 100*np.arange(10)
ED_elec_reqs = 0.5*np.arange(10)

if gather_data:

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
fig.colorbar(pcm1, ax=axs[0], label="Li Recovered (kg/hr)")
axs[0].set_title("Li Recovered (kg/hr)")
axs[0].set_ylabel("Li Price")
axs[0].set_xlabel("ED Electricity Requirement")

pcm2 = axs[1].pcolormesh(ED_elec_reqs, Li_prices, MED, shading='auto')
fig.colorbar(pcm2, ax=axs[1], label="MED Capacity (m³/h)")
axs[1].set_title("MED Capacity (m³/h)")
axs[1].set_xlabel("ED Electricity Requirement")

"""
pcm3 = axs[2].pcolormesh(ED_elec_reqs, Li_prices, RO, shading='auto')
fig.colorbar(pcm3, ax=axs[2], label="RO Capacity (m³/h)")
axs[2].set_title("RO Capacity (m³/h)")
axs[2].set_xlabel("ED Electricity Requirement")
"""
plt.tight_layout()
plt.show()
plt.savefig("physor_updated.png",dpi=1000)

