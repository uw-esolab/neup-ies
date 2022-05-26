import pandas as pd
import scipy as sp
import os
import matplotlib.pyplot as plt 
import numpy as np
import pvlib 
from math import log

# Load all weather files in the directory
fnames = list(filter(lambda x: 'csv' in x, os.listdir(os.getcwd())))

# Compile all years of data into a single dataframe
alldfs = []
for fname in fnames:
    alldfs.append(pd.read_csv(fname, skiprows=2))
df = pd.concat(alldfs, axis=0)

# Report basic stuff
print(df.columns)
print(len(df))


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
csky = tus.get_clearsky(times, model='ineichen', linke_turbidity=1.5).dni
# Add the pvlib computed clear sky DNI to the dataframe in a new column
df['pvlibcs'] = csky.values
# df.plot(y=['pvlibcs','Clearsky DNI'])

# Normalize the observed DNI to clear sky, and compute as deviation from clearsky. Clip values to the range 0..1. Get rid of any nan's.
dnis = (((df["pvlibcs"]-df["DNI"])/df["pvlibcs"]).fillna(0.).clip(0,1))
# Store the scaled DNI value in the dataframe. 
df['DNIs'] = dnis
# dnis.plot()
# df.plot(y=['pvlibcs', 'DNI'])
# plt.show()

# Fourier transform
yf = np.abs(sp.fft.fft(df.DNIs.values))
# Transform the x values
xf = 1./sp.fft.fftfreq(len(df.DNIs), d=0.5)
plt.loglog(xf[1:], yf[1:])
# Come up with some cutoff for significant frequencies from the power density plot. 
cutoff = lambda x: 285.39*log(max(x,1e-17)) + 290.48
plt.loglog(xf[1:], [cutoff(xx) for xx in xf[1:]], 'r--')
plt.show()

# -----------------------------------------------
# From here on, I have no idea what I'm doing.... 


for i in range(len(yf)):
    if yf[i] < cutoff(xf[i]):
        yf[i] = 0

d = np.multiply( sp.fft.ifft(yf) , csky.values )
t = np.arange(0,len(d)/2,0.5)

plt.plot(t[:8760*2],d[:8760*2])
plt.show()