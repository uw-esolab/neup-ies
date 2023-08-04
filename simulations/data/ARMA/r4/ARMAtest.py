# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:47:01 2022

@author: aidan
"""

import pandas as pd
import statsmodels as sm
import numpy as np
import pvlib
import matplotlib.pyplot as plt

"Load All weather data"

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
solar_data = df["DNI"]

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

y = np.zeros(len(csky.values))


#set path to weather file 
filename_orig = 'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r3/solar_data_without_clearsky.csv'

#specify headers of interest
headers = ["Time", "DNI"]

#read csv data into datafram
original_data = pd.read_csv(filename_orig ,names=headers)

#extract data to lists
x = original_data["Time"]
oy = original_data["DNI"]
xa = np.array(x[1:])
xb = xa.astype(float)
oya = np.array(oy[1:])
oyb = oya.astype(float)


arma = sm.arima.model.ARIMA(oyb, order =(1,2,3))
results = arma.fit()
prediction = results.predict(0,len(csky.values))
print(prediction)

# Convert back to DNI. 
for i in range(len(csky.values) -1):
    # Basic transform
    dd = float(csky.values[i]) - float(prediction[i+1])
    # Constrain values to be physical
    dd = max(0,dd)
    dd = min(dd, csky.values[i])
    # Assign result
    y[i] = dd

xa = np.array(x[1:])
xb = xa.astype(float)
ya = np.array(y)
yb = ya.astype(float)

plt.plot(prediction, color='r')
#plt.xlim(0, 25000)
plt.xticks(np.arange(0, 600, 100))
plt.yticks(np.arange(0, 1100, 100))
plt.title("Prediction")
plt.xlabel("Time / minutes")
plt.ylabel("DNI / W/m^2")
plt.show()

plt.plot(xb,yb, color='r')
#plt.xlim(0, 25000)
plt.xticks(np.arange(0, 600000, 100000))
plt.yticks(np.arange(0, 1100, 100))
plt.title("DNI over course of year for Las Vegas - syn data (SV = 8)")
plt.xlabel("Time / minutes")
plt.ylabel("DNI / W/m^2")
plt.show()

plt.plot(xb,yb,color='r')
plt.xlim(0, 20160)
plt.title("DNI over course of first two weeks of a year for Las Vegas - syn data (SV = 8)")
plt.xlabel("Time / minutes")
plt.ylabel("DNI / W/m^2")
#plt.xticks(np.arange(0, 600000, 100000))
#plt.yticks(np.arange(0, 1100, 100))
plt.show()

dni_without_zeros = np.ma.masked_equal(yb,0)

plt.hist(dni_without_zeros, bins=40, edgecolor='black')
plt.show()

# sort the data in ascending order
sorted_dni = np.sort(dni_without_zeros.compressed())
# get the cdf values of y
print(len(sorted_dni))
x_list = np.arange(len(sorted_dni)) / float(len(sorted_dni))
# plotting
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.title('CDF for DNI ')
plt.plot(sorted_dni, x_list, marker='o')