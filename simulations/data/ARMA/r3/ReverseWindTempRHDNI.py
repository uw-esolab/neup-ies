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
solar_data = df["DHI"]

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
cskyn = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).dni
cskyh = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).dhi
cskyg = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).ghi
# Add the pvlib computed clear sky DNI to the dataframe in a new column
df['pvlibcs'] = cskyn.values
# df.plot(y=['pvlibcs','Clearsky DNI'])


filename_origW = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/synDataWind_1.csv"
filename_origT = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/synDataT_1.csv"
filename_origRH = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/synDataRH_1.csv"
filename_origS = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/synDataDG_1.csv"
filename_orig = "C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r4/solar_data_without_clearsky15.csv"
headersW = ["Time", "Wind"]
headersT = ["Time", "Temp"]
headersRH = ["Time","Temp","RH"]
headersS = ["Time","DNI","GHI"]
headers = ["Time", "DNI","DHI","GHI", "Temp", "Precip", "RH", "Wind"]

dataT = pd.read_csv(filename_origT ,names=headersT)
dataW = pd.read_csv(filename_origW ,names=headersW)
dataRH = pd.read_csv(filename_origRH ,names=headersRH)
dataS = pd.read_csv(filename_origS ,names=headersS)
original_data = pd.read_csv(filename_orig ,names=headers)

x = dataT["Time"]
T = dataT["Temp"]
W = dataW["Wind"]
RH = dataRH["RH"]
DNI = dataS["DNI"]
GHI = dataS["GHI"]
ox = original_data["Time"]
oT = original_data["Temp"]
oW = original_data["Wind"]
oRH = original_data["RH"]
oDNI = original_data["DNI"]
oGHI = original_data["GHI"]

x = np.array(x[1:])
x = x.astype(float)

T = np.array(T[1:])
T = T.astype(float)

W = np.array(W[1:])
W = abs(W.astype(float))

DNI = np.array(DNI[1:])
DNI = DNI.astype(float)

GHI = np.array(GHI[1:])
GHI = GHI.astype(float)

#create a second list for the corrected data to be added to
S1 = np.zeros(len(cskyn.values))

# Convert back to DNI. 
for i in range(len(cskyn.values) -1):
    # Basic transform
    dd = float(cskyn.values[i]) - float(DNI[i+1])
    # Constrain values to be physical
    dd = max(0,dd)
    dd = min(dd, cskyn.values[i])
    # Assign result
    S1[i] = dd
    
    #create a second list for the corrected data to be added to
S2 = np.zeros(len(cskyg.values))

# Convert back to DNI. 
for i in range(len(cskyg.values) -1):
    # Basic transform
    dd = float(cskyg.values[i]) - float(GHI[i+1])
    # Constrain values to be physical
    dd = max(0,dd)
    dd = min(dd, cskyg.values[i])
    # Assign result
    S2[i] = dd

RH = np.array(RH[1:])
RH = abs(RH.astype(float))
for i in range(len(RH)):
    if RH[i] > 100:
        RH[i] = 100


ox = np.array(ox[1:])
ox = ox.astype(float)

oT = np.array(oT[1:])
oT = oT.astype(float)

oW = np.array(oW[1:])
oW = oW.astype(float)

oRH = np.array(oRH[1:])
oRH = oRH.astype(float)

oDNI = np.array(oDNI[1:])
oDNI = oDNI.astype(float)

oGHI = np.array(oGHI[1:])
oGHI = oGHI.astype(float)

#create a second list for the corrected data to be added to
oS1 = np.zeros(len(cskyn.values))

# Convert back to DNI. 
for i in range(len(cskyn.values) -1):
    # Basic transform
    dd = float(cskyn.values[i]) - float(oDNI[i+1])
    # Constrain values to be physical
    dd = max(0,dd)
    dd = min(dd, cskyn.values[i])
    # Assign result
    oS1[i] = dd
    
#create a second list for the corrected data to be added to
oS2 = np.zeros(len(cskyg.values))

# Convert back to DNI. 
for i in range(len(cskyg.values) -1):
    # Basic transform
    dd = float(cskyg.values[i]) - float(oGHI[i+1])
    # Constrain values to be physical
    dd = max(0,dd)
    dd = min(dd, cskyg.values[i])
    # Assign result
    oS2[i] = dd
    
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


plt.plot(x,T, color='r')
plt.plot(ox,oT, color='b')
plt.xlim(100000, 120000)
plt.title("Temperature over course of two weeks for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Temp / DegC")
plt.show()

plt.plot(x,T, color='r')
plt.plot(ox,oT, color='b')
#plt.xlim(100000, 120000)
plt.title("Temperature over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Temp / DegC")
plt.show()

plt.plot(x,W, color='r')
plt.plot(ox,oW, color='b')
plt.xlim(100000, 120000)
plt.title("Wind over course of two weeks for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Wind Speed / ms^{-1}")
plt.show()

plt.plot(x,W, color='r')
plt.plot(ox,oW, color='b')
#plt.xlim(100000, 120000)
plt.title("Wind over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Wind Speed / ms^{-1}")
plt.show()

plt.plot(x,RH, color='r')
plt.plot(ox,oRH, color='b')
plt.xlim(100000, 120000)
plt.title("Relative Humidity over course of two weeks for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Relative Humidity")
plt.show()

plt.plot(x,RH, color='r')
plt.plot(ox,oRH, color='b')
#plt.xlim(100000, 120000)
plt.title("Relative Humidity over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("Relative Humidity")
plt.show()


plt.plot(x,S1, color='r')
plt.plot(ox,oS1, color='b')
plt.xlim(100000, 120000)
plt.title("DNI over course of two weeks for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("DNI")
plt.show()

plt.plot(x,S1, color='r')
plt.plot(ox,oS1, color='b')
#plt.xlim(100000, 120000)
plt.title("DNI over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("DNI")
plt.show()

plt.plot(x,S2, color='r')
plt.plot(ox,oS2, color='b')
plt.xlim(100000, 120000)
plt.title("GHI over course of two weeks for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("GHI")
plt.show()

plt.plot(ox,solar_data)
plt.plot(x,S2, color='r')
plt.plot(ox,oS2, color='b')
#plt.xlim(100000, 120000)
plt.title("GHI over course of year for Las Vegas")
plt.xlabel("Time / minutes")
plt.ylabel("GHI")
plt.show()

plt.hist(T, bins=20, edgecolor='black')
plt.hist(oT, bins=20, edgecolor='red')
plt.title("Distribution for Temperature")
plt.show()

plt.hist(W, bins=20, edgecolor='black')
plt.hist(oW, bins=20, edgecolor='red')
plt.title("Distribution for Wind")
plt.show()


plt.hist(RH, bins=20, edgecolor='black')
plt.hist(oRH, bins=20, edgecolor='red')
plt.title("Distribution for RH")
plt.show()

plt.hist(S1, bins=20, edgecolor='black')
plt.hist(oS1, bins=20, edgecolor='red')
plt.title("Distribution for DNI")
plt.show()

plt.hist(oS2, bins=20, edgecolor='red')
plt.hist(S2, bins=20, edgecolor='black')
plt.title("Distribution for GHI")
plt.show()