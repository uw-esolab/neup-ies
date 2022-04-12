import pandas as pd
from scipy.fft import fft, ifft, fftfreq, rfft, rfftfreq, irfft
import os
import matplotlib.pyplot as plt 
import numpy as np
import pvlib 
from math import log
import random

# Load all weather files in the directory
fnames = list(filter(lambda x: 'csv' in x, os.listdir(os.getcwd())))

# Compile all years of data into a single dataframe
alldfs = []
for fname in fnames: #[0:5]:
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
csky = tus.get_clearsky(times, model='ineichen', linke_turbidity=2.5).dni
# Add the pvlib computed clear sky DNI to the dataframe in a new column
df['pvlibcs'] = csky.values
# df.plot(y=['pvlibcs','Clearsky DNI'])

# Do transform on absolute deviation of DNI from clear sky, rather than DNI itself
dnis = (df["pvlibcs"]-df["DNI"]).fillna(0.)
# Store the scaled DNI value in the dataframe. 

print(type(dnis))

a = []

#add noise
for i in dnis:
    i = a.append(i*random.random())


#df['DNIs'] = dnis
df['DNIs'] = a

t = np.arange(0,len(a)/2,0.5)

plt.plot(t[:8760*2],a[:8760*2], label="Synthetic DNI")
#plt.plot(t[:8760*2],df[df.Year == year].DNI, 'r-', label="Actual DNI")
plt.xlim(0, 500)
plt.legend()
plt.show()
#df['DNIs'] = a

# Fourier transform
yf = rfft(df.DNIs.values)
# Transform the x values
xf = rfftfreq(len(df.DNIs), d=0.5)
plt.loglog(xf, np.abs(yf))
# Come up with some cutoff for significant frequencies from the power density plot. 
cutoff = lambda x: 0 #lambda x: 285.39*log(max(x,1e-17)) + 290.48
plt.loglog(xf, [cutoff(xx) for xx in xf], 'r--')
plt.show()

# -----------------------------------------------
# From here on, I have no idea what I'm doing.... 

#add some noise to signal
random.random()*50


for i in range(len(yf)):
    if np.abs(yf[i]) < cutoff(xf[i]):
        yf[i] = 0

# Do the inverse transform
d= np.zeros(len(csky.values))
# raw inverse transform
iyf = irfft(yf)
# Convert back to DNI. 
for i in range(len(csky.values)):
    # Basic transform
    dd = csky.values[i] - iyf[i]
    # Constrain values to be physical
    dd = max(0,dd)
    dd = min(dd, csky.values[i])
    # Assign result
    d[i] = dd 
# Store the synthetic data alongside the other weather data
df['DNIsynth'] = d
# Print some stats
for year in set(df.Year.values):
    print(year)
    print(df[df.Year == year][["DNI","DNIs","DNIsynth"]].describe())

t = np.arange(0,len(d)/2,0.5)

# Plot the data
plt.plot(t[:8760*2],d[:8760*2], label="Synthetic DNI")
plt.plot(t[:8760*2],df[df.Year == year].DNI, 'r-', label="Actual DNI")
plt.legend()
plt.show()

plt.plot(t[:8760*2],d[:8760*2], label="Synthetic DNI")
plt.plot(t[:8760*2],df[df.Year == year].DNI, 'r-', label="Actual DNI")
plt.xlim(0, 500)
plt.legend()
plt.show()