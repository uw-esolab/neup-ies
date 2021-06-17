# Import packages
import pandas as pd
import datetime

# Open  data file
data_file_name = "862903_43.09_-89.42_2018.csv"
data_file = pd.read_csv(data_file_name, header=2).astype(int)

# Extract DHI solar data
target_data = data_file["Temperature"].astype(float)

print(target_data.loc[[0]])

# Extract time data
date=pd.to_datetime(data_file[['Year', 'Month', 'Day', 'Hour', 'Minute']])

# Convert time data to minutes
delta=(pd.to_datetime(data_file['Year'].astype(str)+'-01-01')-date)

time_data=delta.dt.total_seconds().abs() //60

# Concatenate time data and solar data
full_data = pd.concat([time_data, target_data], axis=1, ignore_index=True)

# Write to a new .csv file
full_data.to_csv(r'/home/una/github/postdoc_uw_madison/NE2/sam_dev/neup-ies/simulations/data/ARMA/parsed_data.csv', header=["Time", "Temperature"], index = False)
