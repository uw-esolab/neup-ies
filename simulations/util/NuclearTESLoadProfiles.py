#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:46:53 2022

@author: gabrielsoto
"""

from pylab import rc
import matplotlib.pyplot as plt
import numpy as np
import os, copy, pickle
rc('axes', linewidth=2)
rc('font', weight='bold', size=12)


class NuclearTESLoadProfiles(object):
    """
    The LoadProfileFigures class is SEPARATE from the PostProcessing family of classes.
    This might in the future best be a part of the PostProcessing family... but doing
    this semi-quickly for the Model 1 Energies paper. 
    
    The class sets up all helper methods used to generate overlaid violin plots
    for load profiles of a given run. The given run is specified by:
        - Pref
        - tshours
        - tariff schedule
        - time horizons (both for SSC and Pyomo)
        
    Data from these runs is assumed to be stored in a pickle file and contains 
    data for a full year's worth of simulation. I keep changing file names, so
    the filename is assumed to be an input (a separate script that creates this
    class will handle filename creation). 
    """
    
    def __init__(self, filepath):
        """ Initializes the NuclearTESLoadProfiles module
        """
        self.filepath = filepath
        
        # set some parameters
        self.hrs_per_day      = 24
        self.days_per_week    = 7
        self.weeks_till_a_summer_day = 26
        
        self.summer_start_day = 150
        self.summer_end_day   = 273
        self.year_end_day     = 365

        # some basic arrays
        self.hrs_in_day = np.arange(0, self.hrs_per_day, 1)
        self.Mondays = np.arange(0, self.year_end_day, 7)
        self.Fridays = np.arange(5, self.year_end_day, 7)
        self.Sundays = np.arange(7, self.year_end_day, 7)

        
        # make sure path exists
        assert os.path.exists(self.filepath), "File not found at {0}.".format( self.filepath )
        
        # extracting data as a dictionary from given filepath
        with open(self.filepath, 'rb') as f:
            self.Storage = pickle.load(f)
            
        self.extract_data()
        self.restructure_array_by_hour()
        
        for key in self.Storage.keys():
            array_byHour  = getattr(self, key + "_array_byHour")
            dict_bySeason = self.restructure_dict_by_season_and_weekday( array_byHour )
            setattr(self, key + "_dict", dict_bySeason)
            

    def extract_data(self):
        """ Extracts raw data from Storage 
        
        This method extracts full year's worth of data from storage dictionary
        and saves it as a member attribute of current object.
        """
        
        for key in self.Storage.keys():
            
            # extract key from dictionary
            tmp_array = self.Storage[key]
            
            # check if data has an associated pint unit, remove if so
            if hasattr(tmp_array, 'u'):
                tmp_array = tmp_array.to('GWh') if '_hour' in str(tmp_array.u) else tmp_array

            tmp_array = tmp_array.m if hasattr(tmp_array, 'm') else tmp_array
            
            # save data as object attribute
            setattr(self, key + "_array", tmp_array)
            
            # delete temporary data structure
            del tmp_array
    
    
    def restructure_array_by_hour(self):
        """ Restructures data arrays by hour of day
        
        This method takes data (saved as member arrays of current object)
        and restructures it by hour of the day. That is, for each hour of the day,
        it populates a new array (24 x 365) with values for all, for example,
        midnights throughout the year. 
        """
        
        empty_list = []
        
        for h in self.hrs_in_day:
            
            # creating time slice for all instances of hour h for days in sim
            start_slice = int( h )
            end_slice   = int( self.year_end_day * self.hrs_per_day )
            t_slice     = slice( start_slice, end_slice, 24 )
            
            # slice data
            for key in self.Storage.keys():
                
                if h == 0:
                    setattr(self, key + "_array_byHour", copy.deepcopy(empty_list) )
                
                array_byHour = getattr(self, key + "_array_byHour")
                array_byHour.append( getattr(self, key + "_array")[t_slice].tolist() )
            
        # restructuring as arrays
        for key in self.Storage.keys():
            array_byHour = getattr(self, key + "_array_byHour")
            array_byHour = np.array(array_byHour)
            setattr(self, key + "_array_byHour", array_byHour)
    
    
    def restructure_dict_by_season_and_weekday( self, data_ ):
        """ Restructures existing data arrays into a data dictionary
        
        Data dictionary is split into Winter and Summer data. Those two
        are then split into WeekDay and WeekEnd data. 

        Inputs:
            data_ : ndarry of size 24 x 365

        Outputs:
            data_arrays : dict 

        """
        data_arrays = {}
        
        data_arrays['Winter'] = {}
        data_arrays['Winter']['WeekDay'] = []
        data_arrays['Winter']['WeekEnd'] = []
        
        data_arrays['Summer'] = {}
        data_arrays['Summer']['WeekDay'] = []
        data_arrays['Summer']['WeekEnd'] = []
    
    
        summer_start = self.summer_start_day
        summer_end   = self.summer_end_day
        
        
        for h in self.hrs_in_day:
            # =========================
            # Weekdays
            tmp_winter_WD    = []
            tmp_summer_WD    = []
            
            for m,f in zip(self.Mondays, self.Fridays):
                # slice of data for all hours h -> only grabbing M-F values (as array)
                
                # =============================================================================
                ### Winter  
        
                # considering the first winter months
                if m < summer_start:
                    
                    # checking if summer started mid-week
                    if f <= summer_start:
                        winter_weekday    = data_[h, slice(m,f,1)]
                    else:
                        winter_weekday    = data_[h, slice(m,summer_start,1)]
                        
                    # storing all data for hour h, M-F as a list in list
                    tmp_winter_WD.append( winter_weekday.tolist() )
                
                # considering the second winter months
                elif m > summer_end:
                
                    winter_weekday    = data_[h, slice(m,f,1)]
                    
                    # storing all data for hour h, M-F as a list in list
                    tmp_winter_WD.append( winter_weekday.tolist() )
                
                # =============================================================================
                ### Summer  
                
                else:
                    # checking if summer started mid-week
                    if f <= summer_end:
                        summer_weekday    = data_[h, slice(m,f,1)]
                    else:
                        summer_weekday    = data_[h, slice(m,summer_end,1)]
                        
                    # storing all data for hour h, M-F as a list in list
                    tmp_summer_WD.append( summer_weekday.tolist() )
        
                    
                
            winter_data    = sum(tmp_winter_WD, [])
            data_arrays['Winter']['WeekDay'].append( winter_data )
            
            summer_data    = sum(tmp_summer_WD, [])
            data_arrays['Summer']['WeekDay'].append( summer_data )
            
            # =========================
            # Weekends
            tmp_winter_WE    = []
            tmp_summer_WE    = []
            
            for f,s in zip(self.Fridays, self.Sundays):
                # =============================================================================
                ### Winter  
        
                # considering the first winter months
                if f < summer_start:
                    
                    # checking if summer started mid-week
                    if s <= summer_start:
                        winter_weekend    = data_[h, slice(f,s,1)]
                    else:
                        winter_weekend    = data_[h, slice(f,summer_start,1)]
                        
                    # storing all data for hour h, F-Su as a list in list
                    tmp_winter_WE.append( winter_weekend.tolist() )
                    
                # considering the second winter months
                elif f > summer_end:
                
                    winter_weekend    = data_[h, slice(f,s,1)]
                    
                    # storing all data for hour h, M-F as a list in list
                    tmp_winter_WE.append( winter_weekend.tolist() )
                
                # =============================================================================
                ### Summer  
                
                else:
                    # checking if summer started mid-week
                    if s <= summer_end:
                        summer_weekend    = data_[h, slice(f,s,1)]
                    else:
                        summer_weekend    = data_[h, slice(f,summer_end,1)]
                        
                    # storing all data for hour h, M-F as a list in list
                    tmp_summer_WE.append( summer_weekend.tolist() )
                    
            winter_WEdata    = sum(tmp_winter_WE, [])
            data_arrays['Winter']['WeekEnd'].append( winter_WEdata )
            
            summer_WEdata    = sum(tmp_summer_WE, [])
            data_arrays['Summer']['WeekEnd'].append( summer_WEdata )

        
        
        for season in ["Winter", "Summer"]:
            for dayOfWeek in ["WeekDay", "WeekEnd"]:
                data_arrays[season][dayOfWeek] = np.array( data_arrays[season][dayOfWeek])
                
        return data_arrays
    

    def create_violin_plot(self, ax, data, full_data, v_color='C0', hide_x=False):
        """ Creates a violin plot for each hour of the day
        
        Method repurposed from https://matplotlib.org/stable/gallery/statistics/customized_violin.html
        Creates a violin plot from a list of lists for each of those sub-lists. 
        In our case, the sub-lists are values for a given hour of the day across 
        the full year (e.g., hour 0 from January to December). Also plots
        medians as a white dot, 1st-3rd quartiles in a black rectangle, and the
        full distribution in color specified by input. 
        
        Inputs:
            ax : Pyplot axis object
            data : list of lists 
            v_color : string
        """
        # full violin plot
        vp = ax.violinplot( data, positions=self.hrs_in_day, 
                            widths=0.75, showmeans=False, showmedians=False,
                            showextrema=False)
            
        for pc in vp['bodies']:
            pc.set_facecolor(v_color)
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        
        quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
        whiskers = np.array([
            self.violin_adjacent_values(sorted_array, q1, q3)
            for sorted_array, q1, q3 in zip(np.sort(data,axis=1), quartile1, quartile3)])
        whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
        
        inds = np.arange(0, len(medians) )
        ax.scatter(inds, medians, marker='o', color='black', s=35, zorder=3)
        ax.scatter(inds, medians, marker='o', color='white', s=25, zorder=3)
        ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
        ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
        
        self.modify_violin_axes( ax, data, full_data, hide_x )
        # ax.grid(True)


    def violin_adjacent_values(self, vals, q1, q3):
        """ Helper method for violin plots to clip boundaries
        
        Method taken from https://matplotlib.org/stable/gallery/statistics/customized_violin.html
        """
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
        
        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value


    def modify_violin_axes(self, ax, data, full_data, hide_x=False ):
        
        if not hide_x:
            ax.set_xlabel("Hour of the Day", fontsize=16, fontweight='bold')
            
        ax.set_xlim([-0.5, 23.5])
        ax.set_xticks(self.hrs_in_day)
        
        max_data = []
        min_data = []
        for s in full_data.keys():
            season_data = full_data[s]
            for k in season_data.keys():
                max_data.append(season_data[k].max())
                min_data.append(season_data[k].min())
        
        max_ = np.max(max_data)
        min_ = np.min(min_data)
        upper_lim = max_ * 1.05
        lower_lim = min_ * 0.05
        lower_lim = -1*lower_lim if lower_lim < 0 else lower_lim
        
        n10 = np.round( np.log10( upper_lim ) )
        upper_rounded = np.round( upper_lim / 10**n10) * 10**n10  if n10 > 0 else np.round(upper_lim*10) / 10
        steps = upper_rounded/4 if n10 > 0 else upper_rounded / 3
        tick_range = np.arange(0, upper_rounded, steps)
        
        ax.set_ylim([lower_lim, upper_lim])
        ax.set_yticks( tick_range.tolist() )
        
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)


    def create_tariff_overlay(self, ax, is_winter=True, is_weekday=True, s_color='C0'):
        
        hrs_per_day   = self.hrs_per_day
        days_per_week = self.days_per_week
        weeks_till_summer = self.weeks_till_a_summer_day
        year_end = self.year_end_day
    
        
        M_or_Su = 0 if is_weekday else 6
        neutral_slice = slice(0, hrs_per_day, 1)
        
        lw = 3
        if 'tariffx' in self.filepath:
            
            p_time = self.time_array
            price  = self.price_array

            # indexing the winter default tariff schedule
            
            winter_start_slice  = 0 + M_or_Su * hrs_per_day
            winter_end_slice    = int(winter_start_slice + hrs_per_day)
            winter_slice = slice(winter_start_slice, winter_end_slice, 1)
            
            # indexing the summer default tariff schedule
            summer_start_slice = int( days_per_week * hrs_per_day * weeks_till_summer + M_or_Su * hrs_per_day)
            summer_end_slice   = int( summer_start_slice + hrs_per_day)
            summer_slice = slice(summer_start_slice, summer_end_slice, 1)
            
            ax2 = ax.twinx()
            
            if is_winter:
                ax2.plot(p_time[neutral_slice]-1, price[winter_slice], linewidth= lw, color=s_color, label="Winter Tariff")
            else:
                ax2.plot(p_time[neutral_slice]-1, price[summer_slice], linewidth= lw, color=s_color, label="Summer Tariff")
                ax2.set_ylabel("Price \nMultiplier", fontsize=14, fontweight='bold')
            
            ax2.set_ylim([0.2, 3.5])
            ax2.yaxis.label.set_color(s_color)
            ax2.tick_params(axis='y', colors=s_color, labelsize=14)
            ax2.spines["right"].set_edgecolor(s_color)
            
            
            
        else:

            p_time = self.time_array_byHour
            
            if is_winter:
                price  = np.hstack([ self.price_dict['Winter'][w] for w in ['WeekDay', 'WeekEnd'] ])
                ax.set_ylabel("Price \nMultiplier", labelpad=20, fontsize=12, fontweight='bold')
            else:
                price  = np.hstack([ self.price_dict['Summer'][w] for w in ['WeekDay', 'WeekEnd'] ])
            
            price_mins = price.min(axis=1)
            price_maxs = price.max(axis=1)
            
            price_mean = np.mean(price, axis=1)
            price_std  = np.std( price, axis=1)
            
            price_std_plus  = price_mean + price_std
            price_std_minus = price_mean - price_std
            
            ax.fill_between( self.hrs_in_day, price_mins, price_maxs, color='C0', alpha=0.4 , label='Price Range over 1 Year')
            ax.fill_between( self.hrs_in_day, price_std_minus, price_std_plus, color='C0', alpha=0.6 , label=r'$\pm 1\sigma$')
            ax.plot( self.hrs_in_day, price_mean, 'darkblue', label='Mean Price')
            
            # zeroline
            xzero = np.linspace(-0.5, 23.5, 1000)
            yzero = np.zeros( len(xzero) )
            ax.plot( xzero, yzero, 'k--')
            
            
            # for d in range(len(price.T)):
            #     ax.plot(self.hrs_in_day, price[:,d], color=s_color, alpha=0.5, linewidth=0.5)

            
            ax.set_xlim([-0.5, 23.5])
            ax.set_xticks(self.hrs_in_day)
            
            if is_winter:
                ax.set_ylim([-3, 8.5])
                ax.set_yticks([-2,0,2,4,6,8])
                ax.legend(loc='best', fontsize=12, bbox_to_anchor=(1.17, 1.2),
                          fancybox=True, shadow=True)
                
                peak_x  = 18
                peak_dx = 2.0
                peak_text = r"$-$22.3"
            else:
                ax.set_ylim([-3, 5])
                ax.set_yticks([-2,0,2,4])
                
                peak_x = 21
                peak_dx = 1.55
                peak_text = r"$-$8.5"
                
            
            ax.annotate("", xy=(peak_x, -3), xytext=(peak_x, 0),
                    arrowprops=dict(arrowstyle="->"))
            ax.text(peak_x-peak_dx, -2, peak_text, fontfamily='monospace', fontweight='bold',fontsize=12)
                
