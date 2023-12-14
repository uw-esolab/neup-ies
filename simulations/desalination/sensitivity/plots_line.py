#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 23:46:12 2023

@author: elizabethkeith
"""

import seaborn as sns
import pandas as pd

names_line = []
for i in range(364):
    names_line.append(i)


names_viol = []
for i in range(23):
    names_viol.append(i)
    

dat_line = pd.read_csv('htTES_line2.csv', header=None, names=names_line)
dat_viol = pd.read_csv('htTES_violin2.csv', header=None, names=names_viol)


#sns.lineplot(data=dat_line)
sns.violinplot(data=dat_viol)