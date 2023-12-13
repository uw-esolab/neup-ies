# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 09:12:30 2023

@author: aidan
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from peakdetect import peakdetect
import os
import pvlib
import metpy.calc
from metpy.units import units
import csv
import random as rnd

def weather_data(fnames):
    
    # Compile all years of data into a single dataframe
    alldfs = []
    for fname in fnames: #[0:5]:
        alldfs.append(pd.read_csv(fname, skiprows=2))
    df = pd.concat(alldfs, axis=0)
    
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
    
    return(times)    

def extract_data(fname,file_dict):
    
    headers = list(pd.read_csv(fname, nrows=0))

    data = pd.read_csv(fname ,names=headers)
    
    for n in range(len(headers)):
        
        y = data[headers[n]]
        y = np.array(y[1:])
        y = y.astype(float)
        
        file_dict[headers[n]] = y
    
    
    return(file_dict)

def write_data(file_dict,fname, exfile):
    weatherdf = []
    weatherdf.append(pd.read_csv(exfile, skiprows=2))
    times = pd.concat(weatherdf, axis=0)
    
    with open(fname, 'w', newline='') as myfile:
         wr = csv.writer(myfile,delimiter=',', quoting=csv.QUOTE_MINIMAL)
         wr.writerow(['Source','Location ID','City','State','Country','Latitude','Longitude','Time Zone','Elevation','Local Time Zone','Clearsky DHI Units','Clearsky DNI Units','Clearsky GHI Units','Dew Point Units','DHI Units','DNI Units','GHI Units','Solar Zenith Angle Units','Temperature Units','Pressure Units','Relative Humidity Units','Precipitable Water Units','Wind Direction Units','Wind Speed Units','Cloud Type -15','Cloud Type 0','Cloud Type 1','Cloud Type 2','Cloud Type 3','Cloud Type 4','Cloud Type 5','Cloud Type 6','Cloud Type 7','Cloud Type 8','Cloud Type 9','Cloud Type 10','Cloud Type 11','Cloud Type 12','Fill Flag 0','Fill Flag 1','Fill Flag 2','Fill Flag 3','Fill Flag 4','Fill Flag 5','Surface Albedo Units','Version'])
         wr.writerow(['NSRDB', '102574', '-', '-.1', '-.2', '35.93', '-115.26', '-8', '1021', '-8.1', 'w/m2', 'w/m2.1', 'w/m2.2', 'c', 'w/m2.3', 'w/m2.4', 'w/m2.5', 'Degree', 'c.1', 'mbar', '%', 'cm', 'Degrees', 'm/s', 'N/A', 'Clear', 'Probably Clear', 'Fog', 'Water', 'Super-Cooled Water', 'Mixed', 'Opaque Ice', 'Cirrus', 'Overlapping', 'Overshooting', 'Unknown', 'Dust', 'Smoke', 'N/A.1', 'Missing Image', 'Low Irradiance', 'Exceeds Clearsky', 'Missing CLoud Properties', 'Rayleigh Violation', 'N/A.2', '3.0.6'])
         wr.writerow(['Year', 'Month', 'Day', 'Hour', 'Minute', 'DHI', 'DNI', 'GHI', 'Clearsky DHI', 'Clearsky DNI', 'Clearsky GHI', 'Cloud Type', 'Dew Point', 'Solar Zenith Angle', 'Surface Albedo', 'Wind Speed', 'Precipitable Water', 'Relative Humidity', 'Temperature', 'Pressure'])
         for n in range(len(file_dict["Time"])):
             row_n = [2020,times["Month"][n],times["Day"][n],times["Hour"][n],times["Minute"][n],file_dict["DHI"][n],file_dict["DNI"][n],file_dict["GHI"][n],file_dict["Clearsky DHI"][n],file_dict["Clearsky DNI"][n],file_dict["Clearsky GHI"][n],file_dict["Cloud Type"][n],file_dict["DP"][n],file_dict["Solar Zenith Angle"][n], file_dict["Surface Albedo"][n],file_dict["Wind"][n],file_dict["Precip"][n],file_dict["RH"][n],file_dict["Temp"][n],file_dict["Pressure"][n]]
             wr.writerow(row_n)


def syn_data_correction(file_dict,cskyg,cskyn,cskyh,sz):
    
    Temp = file_dict["Temp"]
    Temp = np.convolve(Temp, np.ones(3)/3, mode='same')
    
    
    # Subtract Cloud Noise from clearsky GHI. 
    GHI = np.zeros(len(cskyg.values))
    for i in range(len(cskyg.values) -1):
        # Basic transform
        GHIi = float(cskyg.values[i]) - file_dict['GHI'][i+1]
        # Constrain values to be physical
        GHIi = max(0,GHIi)
        GHIi = min(GHIi, cskyg.values[i])
        # Assign result
        GHI[i] = GHIi
    GHI = np.convolve(GHI, np.ones(2)/2, mode='same')
    
    DNI = np.zeros(len(cskyn.values))
    # Subtract Cloud Noise from clearsky DNI. 
    for i in range(len(cskyn.values) -1):
        # Basic transform
        DNIi = float(cskyn.values[i]) - file_dict['DNI'][i+1]*1.0
        # Constrain values to be physical
        DNIi = max(0,DNIi)
        DNIi = min(DNIi, cskyn.values[i])
        # Assign result
        DNI[i] = DNIi
    DNI = np.convolve(DNI, np.ones(2)/2, mode='same')
    
    #create empty array for DHI
    DHI = np.zeros(len(cskyh.values))
    #create Envelope for physical DHI
    dhienv = np.zeros(len(cskyn.values))
    max_time = max(file_dict['Time'])
    for i in range(len(cskyh.values) -1):
        dhienv[i] = 400 + 150*np.cos(2*np.pi*((file_dict['Time'][i]-(max_time/2)+25500)/max_time))
    
    # Generate DHI from GHI and DNI 
    for i in range(len(cskyh.values)-1):
        # Basic transform
        #DHIi = cskyh.values[i] - file_dict['GHI'][i+1] + (file_dict['DNI'][i+1])*np.cos(np.radians(sz.values[i]))
        DHIi = GHI[i] - DNI[i]*np.cos(np.radians(sz.values[i]))
        # Constrain values to be physical
        DHIi = max(0,DHIi)
        if cskyh.values[i] == 0:
            DHIi = 0
        # Assign result
        DHI[i] = DHIi
    
    #constric delta DHI and use envelope
    for i in range(len(DHI)-1):
        # delDHI = DHI[i+1]-DHI[i]
        # if delDHI > 150:
        #     DHI[i+1] = DHI[i]+delDHI/20
        if DHI[i] > dhienv[i]:
            DHI[i] = dhienv[i]
    #DHI = np.convolve(DHI, np.ones(2)/2, mode='same')
    
    #add limits to relative humidity
    RH = np.zeros(len(file_dict['RH']))
    for i in range(len(file_dict['RH'])):
        RHi = abs(file_dict['RH'][i])
        if RHi > 100:
            RHi = 200 - RHi - 3*rnd.random()
        if RHi < 3:
            RHi = 3 + 2*rnd.random() - 2*rnd.random()
        RH[i] = RHi
    RH = np.convolve(RH, np.ones(3)/3, mode='same')
        
    #add limits to precipitation humidity
    Precip = np.zeros(len(file_dict['Precip']))
    for i in range(len(file_dict['Precip'])):
        Precipi = abs(file_dict['Precip'][i])
        Precip[i] = Precipi
    Precip = np.convolve(Precip, np.ones(5)/5, mode='same')
        
    #add limits to wind humidity
    Wind = np.zeros(len(file_dict['Wind']))
    for i in range(len(file_dict['Wind'])):
        Windi = abs(file_dict['Wind'][i])
        Wind[i] = Windi
    Wind = np.convolve(Wind, np.ones(5)/5, mode='same')
        
    #generate dew point
    DP = np.zeros(len(file_dict['RH']))
    for i in range(len(file_dict['RH'])):
        if file_dict['Temp'][i] > 35 and file_dict['RH'][i] < 4:
            DPi = metpy.calc.dewpoint_from_relative_humidity(35 * units.degC, 4* units.percent).magnitude
        else: 
            DPi = metpy.calc.dewpoint_from_relative_humidity(file_dict['Temp'][i] * units.degC, file_dict['RH'][i]* units.percent).magnitude
        DP[i] = DPi/1.8
    # for i in range(len(DP)-1):
    #     delDP = DP[i+1]-DP[i]
    #     if abs(delDP) > 5:
    #         DP[i+1] = DP[i]+delDP/10
    DP = np.convolve(DP, np.ones(3)/3, mode='same')
    
    Cloud_type = np.zeros(len(file_dict["Time"]))
    Pressure = np.zeros(len(file_dict["Time"]))
    Solar_Zenith = np.zeros(len(file_dict["Time"]))
    Clearsky_DNI = np.zeros(len(file_dict["Time"]))
    Clearsky_DHI = np.zeros(len(file_dict["Time"]))
    Clearsky_GHI = np.zeros(len(file_dict["Time"]))
    Surface = np.zeros(len(file_dict["Time"]))
    
    for i in range(len(file_dict["Time"])):
        Cloud_type[i] = 1
        Clearsky_DNI[i] = cskyn.values[i]
        Clearsky_DHI[i] = cskyh.values[i]
        Clearsky_GHI[i] = cskyg.values[i]
        Solar_Zenith[i] = sz.values[i]
        Surface[i] = 0.208
        Pressure[i] = 900

    file_dict["GHI"] = GHI
    file_dict["DNI"] = DNI
    file_dict["DHI"] = DHI
    file_dict["RH"] = RH
    file_dict["Temp"] = Temp
    file_dict["Precip"] = Precip
    file_dict["Wind"] = Wind
    file_dict["DP"] = DP
    file_dict["Cloud Type"] = Cloud_type
    file_dict["Clearsky DHI"] = Clearsky_DHI
    file_dict["Clearsky DNI"] = Clearsky_DNI
    file_dict["Clearsky GHI"] = Clearsky_GHI
    file_dict["Pressure"] = Pressure
    file_dict["Surface Albedo"] = Surface
    file_dict["Solar Zenith Angle"] = Solar_Zenith
    
    return(file_dict)

def orig_solar_data_correction(file_dict,cskyg,cskyn,cskyh,sz):
    
    # Subtract Cloud Noise from clearsky GHI. 
    GHI = np.zeros(len(cskyg.values))
    for i in range(len(cskyg.values) -1):
        # Basic transform
        GHIi = float(cskyg.values[i]) - file_dict['GHI'][i+1]
        # Constrain values to be physical
        GHIi = max(0,GHIi)
        #dd = min(dd, cskyg.values[i])
        # Assign result
        GHI[i] = GHIi
    
    DNI = np.zeros(len(cskyn.values))
    # Subtract Cloud Noise from clearsky DNI. 
    for i in range(len(cskyn.values) -1):
        # Basic transform
        DNIi = float(cskyn.values[i]) - file_dict['DNI'][i+1]
        # Constrain values to be physical
        DNIi = max(0,DNIi)
        DNIi = min(DNIi, cskyn.values[i])
        # Assign result
        DNI[i] = DNIi
        
        #generate dew point
    DP = np.zeros(len(file_dict['RH']))
    for i in range(len(file_dict['RH'])):
        DPi = metpy.calc.dewpoint_from_relative_humidity(file_dict['Temp'][i] * units.degC, file_dict['RH'][i]* units.percent).magnitude
        DP[i] = DPi

    file_dict["GHI"] = GHI
    file_dict["DNI"] = DNI
    file_dict["DP"] = DP
    
    return(file_dict)

