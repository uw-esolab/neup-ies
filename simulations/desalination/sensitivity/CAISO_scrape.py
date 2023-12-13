import pandas as pd
import gridstatus
caiso = gridstatus.CAISO()


# Information for CAISO electric prices

# Desired time frame
start_date = pd.Timestamp("Jan 1, 2021").normalize()
end_date = pd.Timestamp("Jan 1, 2022").normalize()

# CAISO node
location = ["DIABLOCN_2_N001"]

# Generating CSV from API
electric_price_schedule = caiso.get_lmp(start=start_date, end=end_date, market="DAY_AHEAD_HOURLY", locations=location) # electric price schedule imported CAISO data for Diablo Canyon NPP
electric_price_schedule.to_csv('electric_price_schedule.csv', index=False)


