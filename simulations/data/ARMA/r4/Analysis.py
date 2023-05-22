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
Variables = ['DHI', 'DNI', 'GHI', 'Dew Point', 'Wind Speed', 'Precipitable Water', 'Relative Humidity', 'Temperature']
Units = ['$W/m^2$','$W/m^2$','$W/m^2$','$^\circ$C','$ms^{-1}$','cm','%','$^\circ$C']
#Variables = ['DHI', 'DNI', 'GHI', 'Wind Speed',]
#Units = ['$W/m^2$','$W/m^2$','$W/m^2$','$ms^{-1}$']
for i in range(1,100,1):
    fullfile = 'TestSyn' + str(i) + '.csv' 
    synfiles.append(fullfile)

for i in range(1998,2020,1):
    fullfile = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/Weather_Data/102574_35.93_-115.26_' + str(i) + '.csv' 
    realfiles.append(fullfile)

for i in range(len(Variables)):
    box_plots(realfiles,synfiles,['102574_35.93_-115.26_2004.csv'],Variables[i], Units[i])

#'102574_35.93_-115.26_2004.csv'