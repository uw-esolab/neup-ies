# Import packages
import pandas as pd
import datetime

# Open  data file
data_file_name = "tucson_az_32.116521_-110.933042_psmv3_60_tmy.csv"
data_file = pd.read_csv(data_file_name, header=2).astype(int)

# Extract DHI solar data
target_data = data_file["Temperature"]

# Extract time data
date=pd.to_datetime(data_file[['Year', 'Month', 'Day', 'Hour', 'Minute']])

# Convert time data to minutes
delta=(pd.to_datetime(data_file['Year'].astype(str)+'-01-01')-date)

time_data=delta.dt.total_seconds().abs() //60

# Concatenate time data and solar data
full_data = pd.concat([time_data, target_data], axis=1, ignore_index=True)

# Write to a new .csv file
full_data.to_csv(r'/home/una/github/postdoc_uw_madison/NE2/sam_dev/neup-ies/simulations/data/ARMA/temp_parsed_data.csv', header=["Time", "Temperature"], index = False)
