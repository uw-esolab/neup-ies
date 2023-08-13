# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:53:34 2022

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from peakdetect import peakdetect
import os
import pvlib

# Load all weather files in the directory
fnames = ['102574_35.93_-115.26_2004.csv']

# Compile all years of data into a single dataframe
alldfs = []
for fname in fnames: #[0:5]:
    alldfs.append(pd.read_csv(fname, skiprows=2))
df = pd.concat(alldfs, axis=0)

# Report basic stuff
print(df.columns)
print(len(df))

# Extract DHI solar data
solar_data = df["Temperature"]

# Extract time data
date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])

# Create a range of dates to evaluate clearsky
dst = df.iloc[0]   #first date in the weather files
dend = df.iloc[-1] #last date in the weather files 
# From the header for the las vegas weather files... 
# 35.93,-115.26,-8,1021,-8
# https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
tus = pvlib.location.Location(35.93,-115.26, 'US/Pacific', 1021, 'Las Vegas')
times = pd.date_range(start='{:04d}-{:02d}-{:02d}'.format(int(dst.Year), int(dst.Month), int(dst.Day)),
                      end='{:04d}-{:02d}-{:02d}'.format(int(dend.Year+1), 1, 1),
                      freq = '30min', 
                      tz=tus.tz)
# Get ride of last time. wraps around and not not needed
times = times.delete(-1)
# Filter out any leap year days, which aren't included in the weather files.
times = times[~((times.month == 2) & (times.day > 28))]

# Calculate clear sky. Use Ineichen model with a low turbidity. Choose value to always exceed observed DNI. 1.5 works for this one.
csky = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).dni
# Add the pvlib computed clear sky DNI to the dataframe in a new column
df['pvlibcs'] = csky.values
# df.plot(y=['pvlibcs','Clearsky DNI'])


filename_orig = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r3/synDataWind_1.csv"
filename_orig2 = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r3/solar_data_without_clearsky.csv"
headers = ["Time", "Temp"]
headers2 = ["Time", "DNI","DHI", "Temp", "Precip", "RH", "Wind"]

original_data = pd.read_csv(filename_orig ,names=headers)
original_data2 = pd.read_csv(filename_orig2 ,names=headers2)

x = original_data["Time"]
z = original_data["Temp"]
x2 = original_data2["Time"]
z2 = original_data2["Temp"]

ox = np.array(x[1:])
oxb = ox.astype(float)
oza = np.array(z[1:])
ozb = oza.astype(float)
ox2 = np.array(x2[1:])
oxb2 = ox2.astype(float)
oza2 = np.array(z2[1:])
ozb2 = oza2.astype(float)

# plt.plot(oxb,oyb, color='r')
# #plt.xlim(0, 25000)
# plt.xticks(np.arange(0, 600000, 100000))
# plt.yticks(np.arange(0, 1100, 100))
# plt.title("DNI Noise Synthetic Signal")
# plt.xlabel("Time / minutes")
# plt.ylabel("DNI / W/m^2")
# plt.show()

# plt.plot(oxb,oyb,color='r')
# plt.xlim(0, 20160)
# plt.title("DNI Noise Synthetic Signal")
# plt.xlabel("Time / minutes")
# plt.ylabel("DNI / W/m^2")
# #plt.xticks(np.arange(0, 600000, 100000))
# #plt.yticks(np.arange(0, 1100, 100))
# plt.show()


xa = np.array(x[1:])
xb = xa.astype(float)
za = np.array(z[1:])
zb = za.astype(float)
xa2 = np.array(x2[1:])
xb2 = xa2.astype(float)
za2 = np.array(z2[1:])
zb2 = za2.astype(float)

plt.plot(xb,zb, color='r')
plt.plot(xb2,zb2, color='b')
plt.xlim(100000, 120000)
plt.title("Temperature over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Temp / DegC")
plt.show()

plt.plot(xb,zb, color='r')
plt.plot(xb2,zb2, color='b')
#plt.xlim(100000, 120000)
plt.title("Temperature over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Temp / DegC")
plt.show()


plt.hist(zb, bins=20, edgecolor='black')
plt.hist(zb2, bins=20, edgecolor='red')
plt.show()