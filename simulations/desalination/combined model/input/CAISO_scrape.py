import pandas as pd
import gridstatus
caiso = gridstatus.CAISO()


# specify desired time frame

date_start = pd.Timestamp("Jan 1, 2021").normalize()
date_end   = pd.Timestamp("Jan 1, 2022").normalize()



# specify desired CAISO node 

node_str = 'DIABLOCN_2_N001'
node     = [node_str]



# Generate csv from CAISO API

sched_elec = caiso.get_lmp(start=date_start, end=date_end, market="DAY_AHEAD_HOURLY", locations=node) 
sched_elec.to_csv(node_str + '_PRICES.csv', index=False)

