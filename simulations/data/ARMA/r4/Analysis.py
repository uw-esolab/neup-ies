# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 08:58:59 2023

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pvlib
from scipy import stats
import statistics
from statsmodels.distributions.empirical_distribution import ECDF
from analysis_fuctions import box_plots

synfiles = []
realfiles = []
for i in range(1,16,1):
    fullfile = 'TestSyn' + str(i) + '.csv' 
    synfiles.append(fullfile)

for i in range(1998,2020,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/Weather_Data/102574_35.93_-115.26_' + str(i) + '.csv' 
    realfiles.append(fullfile)

box_plots(realfiles,synfiles,['102574_35.93_-115.26_2004.csv'],'DNI')