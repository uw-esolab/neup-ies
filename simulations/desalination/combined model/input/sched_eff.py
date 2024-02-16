import pandas as pd


# specify CAISO node 

node = 'DIABLOCN_2_N001'


# read csv file with temperature data

sched_temp     = node + '_TEMPS.csv'
sched_temp     = pd.read_csv(sched_temp, usecols=['Temperature'])
sched_temp     = sched_temp['Temperature'].values


# create list of ambient temperature efficiency correction indexed to hour of the year

sched_eff      = []                                                                                                                                 #efficiency schedule as function of ambient temperatures
for temp in sched_temp:
    eff  = (-0.0036*temp + 1.2161) - 0.073
    sched_eff.append(eff)
    

# write csv file with efficiency correction data

df_out    = pd.DataFrame(sched_eff, columns=['Efficiency'])
df_out.to_csv('sched_eff.csv')



