#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 21:43:50 2023

@author: elizabethkeith
"""

import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
import pandas as pd

data1     = pd.read_csv('results1.csv')         
data2     = pd.read_csv('results2.csv')         
data3     = pd.read_csv('results3.csv')         
# data4   = pd.read_csv('results4.csv')         
# data5   = pd.read_csv('results5.csv')         


# electricity prices & power output
fig1, ax1 = plt.subplots()
lns1      = ax1.plot(data3.t, data3.P_elec, linewidth=0.5, color='red', label='Electricity Prices')
ax2       = ax1.twinx()
lns2      = ax2.plot(data3.t, data3.w_dot, linewidth=0.5, color='blue', label='Power Output')
ax1.set_xlabel("Timestep (hr)")
ax1.set_ylabel("Electricity Prices ($/MWh)")
ax2.set_ylabel("Power Output (MW)")
plt.xlim(3200,3300)

lns_a     = lns1 + lns2
labs_a    = [l.get_label() for l in lns_a]
ax1.legend(lns_a, labs_a, loc='upper right')


# extraction mass flow & power output
fig2, ax3 = plt.subplots()
lns3      = ax3.plot(data3.t, data3.m_dot_e, linewidth=0.1, color='red')
ax4       = ax3.twinx()
lns4      = ax4.plot(data3.t, data3.w_dot, linewidth=0.1, color='blue')
ax3.set_xlabel('Timesteps (hr)')
ax3.set_ylabel('Extraction from Steam Cycle (kg/s)')
ax4.set_ylabel('Power Output (MW)')

lns_b     = lns3 + lns4
labs_b    = [l.get_label() for l in lns_b]
ax3.legend(lns_b, labs_b, loc='upper right')



# discharge from cycle TES & power production
fig3, ax5 = plt.subplots()
lns5      = ax5.plot(data3.t, data3.m_dot_cs, linewidth=0.1, color='red')
ax6       = ax5.twinx()
lns6      = ax6.plot(data3.t, data3.w_dot, linewidth=0.1, color='blue')
ax5.set_xlabel('Timesteps (hr)')
ax5.set_ylabel('Mass Flow from Cyle TES (kg/s)')
ax6.set_ylabel('Power Output (MW)')

lns_c     = lns5 + lns6
labs_c    = [l.get_label() for l in lns_c]
ax3.legend(lns_c, labs_c, loc='upper right')


# distillate production & power production
fig4, ax7 = plt.subplots()
lns7      = ax7.plot(data3.t, data3.v_dot, linewidth=0.5, color='red', label='Distillate')
ax8       = ax7.twinx()
lns8      = ax8.plot(data3.t, data3.w_dot, linewidth=0.5, color='blue', label='Electricity')
ax7.set_xlabel('Timesteps (hr)')
ax7.set_ylabel('Distillate Produced (kg/s)')
ax8.set_ylabel('Power Produced (MW)')
plt.xlim(3200,3300)

lns_d     = lns7 + lns8
labs_d    = [l.get_label() for l in lns_d]
ax7.legend(lns_d, labs_d, loc='upper right')



# cycle storage & power production
fig5, ax9 = plt.subplots()
lns9      = ax9.plot(data3.t, data3.m_ch, linewidth=0.5, color='red', label='Mass in High Temp TES')
ax10      = ax9.twinx()
lns10     = ax10.plot(data3.t, data3.w_dot, linewidth=0.5, color='blue', label='Power Produced')
lns11     = ax10.plot(data3.t, data3.P_elec, linewidth=0.5, color='green', label='Electricity Prices')
ax9.set_xlabel('Timesteps (s)')
ax9.set_ylabel('Mass in Cycle Storage (kg)')
ax10.set_ylabel('Power Produced (MW)')
plt.xlim(3200,3300)

lns_e     = lns9 + lns10 + lns11
labs_e    = [l.get_label() for l in lns_e]
ax9.legend(lns_e, labs_e, loc='upper right')







