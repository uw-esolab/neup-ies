
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 600

tit = '../../../runs/base/base.csv'
dat = pd.read_csv(tit)


sched_elec = pd.read_csv('../../../input/DIABLOCN_2_N001_PRICES.csv')
sched_elec = sched_elec['LMP'].values
sched_elec = pd.Series(sched_elec)


t      = dat.t
ctes   = (dat.m_ch*1.5*(565-330)/3600)/1000     # convert mass to energy units
dtes   = dat.m_dh*421.3/3600/1000               # convert mass to energy units
power  = dat.w_dot               
dist   = dat.v_dot*86.4/100                     # convert from kg/s to 100's m^3/day




plt.figure(figsize=(8,12),dpi=1000)


# first subplot: energy in thermal storage
x1 = plt.subplot(3, 1, 1)
lns1 = ctes.plot(color='darkred', label='High Temperature', linewidth=4)
x2 = plt.subplot(3, 1, 1)
lns2 = dtes.plot(color='royalblue', label='Low Temperature', linewidth=4)


# second subplot: power production and distillate production
x3 = plt.subplot(3, 1, 2)
lns3 = power.plot(ax=plt.gca(), color='darkorange', label='Electricity', linewidth=4)
x4 = plt.subplot(3, 1, 2)
lns4 = dist.plot(ax=plt.gca(), color='teal', label='Distillate', linewidth=4)


# third subplot: CAISO electricity prices
x5 = plt.subplot(3, 1, 3)
lns5 = sched_elec.plot(ax=plt.gca(), color='pink', label='Price', linewidth=4)


# assign characteristics to first subplot
x1.set_xticks([])
x1.set_xlim([6216,6288])
x1.set_ylabel('Thermal Storage Inventory [MWh]')
x1.legend(loc='upper right')
x1.set_title('Dispatch Profile: Base Case')


# assign characteristics to second subplot
x3.set_xticks([])
x3.set_xlim([6216,6288])
x3.set_ylabel('Electricity Production [MW]\n dist Production [10^2 m^3/day]')
x3.legend(loc='upper right')


# assign characteristics to third subplot
x5.set_xticks([6216, 6240, 6264, 6288])
x5.set_xlim([6216, 6288])
x5.set_ylim([-10, 110])
x5.set_ylabel('Electricity Prices ($/MWh)')
x5.set_xticklabels(['Sept 16 \n 12:00 AM', 'Sept 17 \n 12:00 AM', 'Sept 18 \n 12:00 AM', 'Sept 19 \n 12:00 AM'])
x5.legend(loc='upper right')



# turn gridlines on
x1.grid(True, axis='both', which='both')
x3.grid(True, axis='both', which='both')
x5.grid(True, axis='y', which='both')






