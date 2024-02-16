#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 12:19:46 2024

@author: elizabethkeith
"""


import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn

plt.figure()
dat_elec = pd.read_csv('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/ltTES_line_elec.csv', names=['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%'], header=None)
elec = sn.lineplot(data=dat_elec)
elec.set_xlim(1500, 1800)
elec.set(xlabel='Timestep (hr)', ylabel='Power Produced (MW)')

plt.figure()
dat_dist = pd.read_csv('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/ltTES_line_dist.csv', names=['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%'], header=None)
dist = sn.lineplot(data=dat_dist)
dist.set_xlim(1500, 1800)
dist.set(xlabel='Timestep (hr)', ylabel='Distillate Produced (kg/s)')

plt.figure()
dat_cstor = pd.read_csv('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/ltTES_line_cstor.csv', names=['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%'], header=None)
dat_cstor = dat_cstor*(1.5*(565-330)/3600)/1000
cstor = sn.lineplot(data=dat_cstor)
cstor.set_xlim(1500, 1800)
cstor.set(xlabel='Timestep (hr)', ylabel='High Temp Inventory (MWh)')


plt.figure()
dat_dstor = pd.read_csv('/Users/elizabethkeith/neup-ies/simulations/desalination/sensitivity/output files/ltTES_line_dstor.csv', names=['25%', '50%', '75%', '100%', '125%', '150%', '175%', '200%'], header=None)
dat_dstor = dat_dstor*421.3/3600/1000
dstor = sn.lineplot(data=dat_dstor)
dstor.set_xlim(1500, 1800)
dstor.set(xlabel='Timestep (hr)', ylabel='Low Temp Inventory (MWh)')
