#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:06:55 2022

@author: gabrielsoto
"""

import os,sys
sys.path.append('..')
import modules.NuclearTES as NuclearTES
from util.FileMethods import FileMethods
import matplotlib.pyplot as plt
import os, pint, time, copy, pickle
import numpy as np
import pyomo.environ as pe
from pylab import rc
rc('axes', linewidth=2)
rc('font', weight='bold',size=12)
u = pint.UnitRegistry()

pid = os.getpid()
print("PID = ", pid)

# =============================================================================
#  Simulation Parameters
# =============================================================================

# modifying inputs
json = "model1_Hamilton_560_tariffx1"   # model1_CAISO_Hamilton # model1_Hamilton_560_tariffx2 # model1_Hamilton_560_tariffx1
dispatch = True
run_loop = True
sscH    = 24   # (hr)
pyoH    = 48   # (hr)
Pref    = 700  # (MW)
tshours = 2    # (hr)

# =============================================================================
#  Extracting Information
# =============================================================================

# locating output directory
output_dir = os.path.join( FileMethods.output_dir, "model1_energies_paper" )
filename = 'loadProfile_PySAM__{0}__2022_02__pyomo_{1:.0f}__horizon_{2:.0f}_{3:.0f}__TES_{4}__PC_{5}.nuctes'.format(
                json, dispatch, sscH, pyoH, tshours, Pref )
NTPath = os.path.join(output_dir, filename)

# pickling
with open(NTPath, 'rb') as f:
    Storage = pickle.load(f)
    
# storage dictionary
gen     = Storage['gen'].m      
q_nuc   = Storage['q_nuc_thermal'].m  
q_pb    = Storage['q_pb'].m      
q_su    = Storage['q_dot_pc_su'].m  
etes    = Storage['e_ch_tes'].m     
price   = Storage['price']   
time    = Storage['time']   

# =============================================================================
#  Helper method to adjust limits of violin plot taken from 
#    https://matplotlib.org/stable/gallery/statistics/customized_violin.html
# =============================================================================

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

# =============================================================================
# Plotting parameters for loops
# =============================================================================
hrs_per_day   = 24

# number of days to loop through
summer_start = 150
summer_end   = 273
year_end = 365

# time arrays, converting to right units and setting arrays to plot
p_time = time.to('hr').m  
hours = np.arange(0, hrs_per_day, 1)
hrs_in_year = np.arange(0, year_end*24, 1)

# slices for the year
winter1 = slice(           0, summer_start, 1)
summer1 = slice(summer_start,   summer_end, 1)
winter2 = slice(  summer_end,     year_end, 1)

# empty arrays to slice through hours of the day
hrs_data   = []
gen_data   = []
q_nuc_data = []
q_pb_data  = []
q_su_data  = []
etes_data  = []

# =============================================================================
# Create Violin Plot
# =============================================================================

def create_violin_plot(ax, data, v_color='C0'):
    
    # full violin plot
    vp = ax.violinplot( data, positions=hours, 
                        widths=0.75, showmeans=False, showmedians=False,
                        showextrema=False)
        
    for pc in vp['bodies']:
        pc.set_facecolor(v_color)
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    
    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(np.sort(data,axis=1), quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    
    inds = np.arange(0, len(medians) )
    ax.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
    
    
def create_tariff_overlay(ax, is_winter=True, is_weekday=True, s_color='C0'):
    
    if 'tariffx' in json:
        days_per_week = 7
        hrs_per_day   = 24
        weeks_till_summer = 26
        
        neutral_slice = slice(0, hrs_per_day, 1)
        
        # indexing the winter default tariff schedule
        M_or_Su = 0 if is_weekday else 6
        winter_start_slice  = 0 + M_or_Su * hrs_per_day
        winter_end_slice    = int(winter_start_slice + hrs_per_day)
        winter_slice = slice(winter_start_slice, winter_end_slice, 1)
        
        # indexing the summer default tariff schedule
        summer_start_slice = int( days_per_week * hrs_per_day * weeks_till_summer + M_or_Su * hrs_per_day)
        summer_end_slice   = int( summer_start_slice + hrs_per_day)
        summer_slice = slice(summer_start_slice, summer_end_slice, 1)
        
        ax2 = ax.twinx()
        
        if is_winter:
            ax2.plot(p_time[neutral_slice]-1, price[winter_slice], linewidth= 4, color=s_color, label="Winter Tariff")
        else:
            ax2.plot(p_time[neutral_slice]-1, price[summer_slice], linewidth= 4, color=s_color, label="Summer Tariff")
            ax2.set_ylabel("Price Multiplier", fontweight='bold')
        
        ax2.set_ylim([0.2, 3.5])
        ax2.yaxis.label.set_color(s_color)
        ax2.tick_params(axis='y', colors=s_color)
        ax2.spines["right"].set_edgecolor(s_color)
        
        
        
    else:
        days_per_week = 7
        hrs_per_day   = 24
        
        neutral_slice = slice(0, hrs_per_day, 1)
        
        M_or_Su = 0 if is_weekday else 6
        
        ax2 = ax.twinx()
        
        for w in range(52):
            
            start_slice = int( days_per_week * hrs_per_day * w)
            end_slice   = int( start_slice + days_per_week * hrs_per_day)
            t_slice     = slice(start_slice, end_slice, 1)
            
            if w == 0:
                s_slice = t_slice
                
            LMP_slice = price[t_slice]
            
            ax2.plot(p_time[neutral_slice], LMP_slice, color='C0', linewidth= 1, alpha=0.5)


# =============================================================================
# Slicing per Hour
# =============================================================================

# loop through each hour of the day [0 -> 24]
#    so e.g. for hour=0, we're looking at every midnight from start of sim
#            til the number of "days" specified above
for h in hours:
    
    #++++++ WINTER - to - SUMMER ========#
    # creating time slice for all instances of hour h for days in sim
    start_slice = int( h )
    end_slice   = int( year_end * 24 )
    t_slice     = slice( start_slice, end_slice, 24 )
    
    # slice data
    hrs_data.append(  hrs_in_year[t_slice] )
    gen_data.append(   gen[t_slice] )
    q_nuc_data.append( q_nuc[t_slice] )
    q_pb_data.append(  q_pb[t_slice] )
    q_su_data.append(  q_su[t_slice] )
    etes_data.append(  etes[t_slice] )
    
# restructuring as arrays
hrs_data   = np.array(hrs_data)
gen_data   = np.array(gen_data)
q_nuc_data = np.array(q_nuc_data)
q_pb_data  = np.array(q_pb_data)
q_su_data  = np.array(q_su_data)
etes_data  = np.array(etes_data)


# =============================================================================
# Weekday vs Weekend Behavior
# =============================================================================

Mondays = np.arange(0, year_end, 7)
Fridays = np.arange(5, year_end, 7)
Sundays = np.arange(7, year_end, 7)


def generate_data_array( data_ ):
    data_arrays = {}
    
     
    data_arrays['Summer'] = {}
    data_arrays['Summer']['WeekDay'] = {}
    data_arrays['Summer']['WeekDay']['data'] = []
    data_arrays['Summer']['WeekDay']['hour'] = []
    data_arrays['Summer']['WeekEnd'] = {}
    data_arrays['Summer']['WeekEnd']['data'] = []
    data_arrays['Summer']['WeekEnd']['hour'] = []
    
    data_arrays['Winter'] = {}
    data_arrays['Winter']['WeekDay'] = {}
    data_arrays['Winter']['WeekDay']['data'] = []
    data_arrays['Winter']['WeekDay']['hour'] = []
    data_arrays['Winter']['WeekEnd'] = {}
    data_arrays['Winter']['WeekEnd']['data'] = []
    data_arrays['Winter']['WeekEnd']['hour'] = []
    
    
    for h in hours:
        # =========================
        # Weekdays
        tmp_winter_WD    = []
        tmp_winter_WD_hr = []
        
        tmp_summer_WD    = []
        tmp_summer_WD_hr = []
        
        for m,f in zip(Mondays, Fridays):
            # slice of data for all hours h -> only grabbing M-F values (as array)
            
            # =============================================================================
            ### Winter  
    
            # considering the first winter months
            if m < summer_start:
                
                # checking if summer started mid-week
                if f <= summer_start:
                    winter_weekday    = data_[h, slice(m,f,1)]
                    winter_hour       = hrs_data[h, slice(m,f,1)]
                else:
                    winter_weekday    = data_[h, slice(m,summer_start,1)]
                    winter_hour       = hrs_data[h, slice(m,summer_start,1)]
                    
                # storing all data for hour h, M-F as a list in list
                tmp_winter_WD.append( winter_weekday.tolist() )
                tmp_winter_WD_hr.append( winter_hour.tolist() )
            
            # considering the second winter months
            elif m > summer_end:
            
                winter_weekday    = data_[h, slice(m,f,1)]
                winter_hour       = hrs_data[h, slice(m,f,1)]
                
                # storing all data for hour h, M-F as a list in list
                tmp_winter_WD.append( winter_weekday.tolist() )
                tmp_winter_WD_hr.append( winter_hour.tolist() )
            
            # =============================================================================
            ### Summer  
            
            else:
                # checking if summer started mid-week
                if f <= summer_end:
                    summer_weekday    = data_[h, slice(m,f,1)]
                    summer_hour       = hrs_data[h, slice(m,f,1)]
                else:
                    summer_weekday    = data_[h, slice(m,summer_end,1)]
                    summer_hour       = hrs_data[h, slice(m,summer_end,1)]
                    
                # storing all data for hour h, M-F as a list in list
                tmp_summer_WD.append( summer_weekday.tolist() )
                tmp_summer_WD_hr.append( summer_hour.tolist() )
    
                
            
        winter_data    = sum(tmp_winter_WD, [])
        winter_hour    = sum(tmp_winter_WD_hr, [])
        data_arrays['Winter']['WeekDay']['data'].append( winter_data )
        data_arrays['Winter']['WeekDay']['hour'].append( winter_hour )
        
        summer_data    = sum(tmp_summer_WD, [])
        summer_hour    = sum(tmp_summer_WD_hr, [])
        data_arrays['Summer']['WeekDay']['data'].append( summer_data )
        data_arrays['Summer']['WeekDay']['hour'].append( summer_hour )
        
        # =========================
        # Weekends
        tmp_winter_WE    = []
        tmp_winter_WE_hr = []
        
        tmp_summer_WE    = []
        tmp_summer_WE_hr = []
        
        for f,s in zip(Fridays, Sundays):
            # =============================================================================
            ### Winter  
    
            # considering the first winter months
            if f < summer_start:
                
                # checking if summer started mid-week
                if s <= summer_start:
                    winter_weekend    = data_[h, slice(f,s,1)]
                    winter_WEhour     = hrs_data[h, slice(f,s,1)]
                else:
                    winter_weekday    = data_[h, slice(f,summer_start,1)]
                    winter_WEhour     = hrs_data[h, slice(f,summer_start,1)]
                    
                # storing all data for hour h, F-Su as a list in list
                tmp_winter_WE.append( winter_weekend.tolist() )
                tmp_winter_WE_hr.append( winter_WEhour.tolist() )
                
            # considering the second winter months
            elif f > summer_end:
            
                winter_weekend    = data_[h, slice(f,s,1)]
                winter_WEhour     = hrs_data[h, slice(f,s,1)]
                
                # storing all data for hour h, M-F as a list in list
                tmp_winter_WE.append( winter_weekend.tolist() )
                tmp_winter_WE_hr.append( winter_WEhour.tolist() )
            
            # =============================================================================
            ### Summer  
            
            else:
                # checking if summer started mid-week
                if s <= summer_end:
                    summer_weekend    = data_[h, slice(f,s,1)]
                    summer_WEhour     = hrs_data[h, slice(f,s,1)]
                else:
                    summer_weekend    = data_[h, slice(f,summer_end,1)]
                    summer_WEhour     = hrs_data[h, slice(f,summer_end,1)]
                    
                # storing all data for hour h, M-F as a list in list
                tmp_summer_WE.append( summer_weekend.tolist() )
                tmp_summer_WE_hr.append( summer_WEhour.tolist() )
                
        winter_WEdata    = sum(tmp_winter_WE, [])
        winter_WEhour    = sum(tmp_winter_WE_hr, [])
        data_arrays['Winter']['WeekEnd']['data'].append( winter_WEdata )
        data_arrays['Winter']['WeekEnd']['hour'].append( winter_WEhour )
        
        summer_WEdata    = sum(tmp_summer_WE, [])
        summer_WEhour    = sum(tmp_summer_WE_hr, [])
        data_arrays['Summer']['WeekEnd']['data'].append( summer_WEdata )
        data_arrays['Summer']['WeekEnd']['hour'].append( summer_WEhour )
    
    
    for season in ["Winter", "Summer"]:
        for dayOfWeek in ["WeekDay", "WeekEnd"]:
            
            data_arrays[season][dayOfWeek]["data"] = np.array( data_arrays[season][dayOfWeek]["data"] )
            data_arrays[season][dayOfWeek]["hour"] = np.array( data_arrays[season][dayOfWeek]["hour"] )
            
    return data_arrays

# =============================================================================
# create figure
# =============================================================================
gen_data_array   = generate_data_array(gen_data)
tch_data_array   = generate_data_array(etes_data)
# q_nuc_data = np.array(q_nuc_data)
# q_pb_data  = np.array(q_pb_data)
# q_su_data  = np.array(q_su_data)
# etes_data  = np.array(etes_data)


def common_plot_things( ax, is_gen=True ):
    ax.set_xlabel("Hour of the Day", fontweight='bold')
    
    ax.set_xlim([-0.5, 23.5])
    ax.set_xticks(hours)
    
    if is_gen:
        ax.set_ylim([-5, Pref*1.1])
    else:
        ax.set_ylim([-5, Pref*tshours*2.2])
    
    
# create full figure

gen_fig = plt.figure(figsize=(18,7))
gen_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
             fontsize=16, fontweight='bold')

tch_fig = plt.figure(figsize=(18,7))
tch_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
             fontsize=16, fontweight='bold')

count = 0
for dayOfWeek in ["WeekDay", "WeekEnd"]:
    for season in ["Winter", "Summer"]:
    
        # defining subplot of full figure (2x2)
        subplot = "22" + str(count+1)
        
        # creating current axis
        gen_ax__ = gen_fig.add_subplot( eval(subplot) )
        tch_ax__ = tch_fig.add_subplot( eval(subplot) )


        # booleans for tariff overlay
        is_winter  = season is "Winter"
        is_weekday = dayOfWeek is "WeekDay"
        
        # create violin plot and tariffs on twined axis
        plotting_gen_array = gen_data_array[season][dayOfWeek]['data'].tolist()
        plotting_tch_array = tch_data_array[season][dayOfWeek]['data'].tolist()
        
        create_violin_plot(    gen_ax__, plotting_gen_array, v_color='C9' )
        create_violin_plot(    tch_ax__, plotting_tch_array, v_color='C4' )
        
        if 'tariffx' in json:
            create_tariff_overlay( gen_ax__, is_winter, is_weekday, s_color='C1')
            create_tariff_overlay( tch_ax__, is_winter, is_weekday, s_color='C1')
        
        # some basic plotting labels
        common_plot_things( gen_ax__ )
        common_plot_things( tch_ax__, is_gen=False )
        
        if count == 0 or count == 1:
            gen_ax__.set_title(season, fontsize=16, fontweight='bold')
            tch_ax__.set_title(season, fontsize=16, fontweight='bold')
        if count == 0 or count == 2:
            gen_ax__.set_ylabel(dayOfWeek + '\n$\\regular_{Power\ Generated \ (MWe)}}$', fontsize=16, fontweight='bold')
            tch_ax__.set_ylabel(dayOfWeek + '\n$\\regular_{Power\ Generated \ (MWe)}}$', fontsize=16, fontweight='bold')
    
        count += 1
        
gen_fig.tight_layout()
tch_fig.tight_layout()


# else:
#     gen_fig = plt.figure(figsize=(18,5))
#     gen_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
#                  fontsize=16, fontweight='bold')
    
#     tch_fig = plt.figure(figsize=(18,5))
#     tch_fig.suptitle("{0} \nDesign Output = {1} MWe, \nTES = {2} hrs".format(json, Pref, tshours), 
#                  fontsize=16, fontweight='bold')
    
#     gen_ax__ = gen_fig.add_subplot( 111 )
#     tch_ax__ = tch_fig.add_subplot( 111 )
    
#     create_violin_plot(    gen_ax__, gen_data.tolist(), v_color='C9' )
#     # create_tariff_overlay( gen_ax__, is_winter, is_weekday, s_color='C1')
    
#     create_violin_plot(    tch_ax__, etes_data.tolist(), v_color='C4' )
#     # create_tariff_overlay( tch_ax__, is_winter, is_weekday, s_color='C1')
    