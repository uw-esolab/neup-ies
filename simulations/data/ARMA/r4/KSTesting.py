# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 09:21:16 2022

@author: aidan
"""

from scipy import stats
import numpy as np

rng = np.random.default_rng()
n1 = 200  # size of first sample
n2 = 300  # size of second sample

#Create two random samples from a normal distribution at loc = loc and size = scale
rvs1 = stats.norm.rvs(size=n1, loc=0., scale=1, random_state=rng)
rvs2 = stats.norm.rvs(size=n2, loc=5, scale=3, random_state=rng)

#perform KS test on two samples
print(stats.ks_2samp(rvs1, rvs2))