# Import packages
import pandas as pd
import datetime
import pvlib 
import matplotlib.pyplot as plt 
import numpy as np

# Open  data file
data_file_name = "862903_43.09_-89.42_2018.csv"
df = pd.read_csv(data_file_name, header=2).astype(int)

# Extract DHI solar data
solar_data = df["DNI"]

# Extract time data
date=pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour', 'Minute']])

# Convert time data to minutes
delta=(pd.to_datetime(df['Year'].astype(str)+'-01-01')-date)

# Create a range of dates to evaluate clearsky
dst = df.iloc[0]   #first date in the weather files
dend = df.iloc[-1] #last date in the weather files 
# From the header for the las vegas weather files... 
# 43.09	-89.42	-6	260	-6
# https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
tus = pvlib.location.Location(43.09,-89.42, 'US/Pacific', 260, 'Madison')
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

time_data=delta.dt.total_seconds().abs() //60

# Concatenate time data and solar data
full_data = pd.concat([time_data, solar_data], axis=1, ignore_index=True)

# Write to a new .csv file
full_data.to_csv(r'C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/solar_data_parsed_DNI.csv', header=["Time", "DNI"], index = False)
# Write solar data with new time points

t = np.arange(0,len(dnis)/2,0.5)

plt.plot(t[:8760*2],dnis[:8760*2], label="Synthetic DNI")
#plt.plot(t[:8760*2],df[df.Year == year].DNI, 'r-', label="Actual DNI")
plt.xlim(0, 500)
plt.legend()
plt.show()