# Convert from minutes into hour angle, omega
# Calculate normalised hour angle, omega_bar, where omega_bar = omega / omega_max for each day
 # Import modules
import pandas as pd
import matplotlib.pyplot as plt
import pvlib
from datetime import datetime

# Read in data from .csv file for time and DNI
input_data = pd.read_csv("/home/una/github/postdoc_uw_madison/NE2/sam_dev/neup-ies/simulations/data/ARMA/solar_data_parsed.csv")
df = pd.DataFrame()
df["Time"]= input_data["Time"]
df["DHI"] = input_data["DHI"]


# Convert to DatetimeIndex for use in pvlib
start_date = pd.Timestamp(2018,1,1).tz_localize('America/Chicago')
df.index = start_date + pd.to_timedelta(df['Time'], unit='min')


# Convert to DatetimeIndex for use in pvlib
#df["Time"] = pd.to_datetime(df["Time"],unit="m",origin="01-01-2018",tz="America/Chicago")
#df.index = pd.DatetimeIndex(df["Time"])

# Calculate equation of time for hour angle calculation - 1st day of the year
equation_of_time = pvlib.solarposition.equation_of_time_spencer71(1)

# Find hour angle
df["Hour angle"] = pvlib.solarposition.hour_angle(df.index, -89.401230, equation_of_time)

df.to_csv(r'/home/una/github/postdoc_uw_madison/NE2/sam_dev/neup-ies/simulations/data/ARMA/hour_angle_parsed_data.csv')
