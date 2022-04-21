 # Import modules
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import Divider, Size

extract_data = False
plot_data = True

 # Set file paths
dirPath = "/home/una/github/postdoc_uw_madison/NE2/sam_dev/neup-ies/simulations/data/ARMA/r1/"
filePrefix = "synData_"
pivotParameter = "Time"
variable = "Temperature"
histories = 100 # Number of RAVEN histories run
outputFile = "arma_output.csv"
original_data = "temperature_data_parsed.csv"
avg_orig_data = "temperature_data_parsed_avg.csv"

if extract_data == True:
    # Extract time data from first file, write to new collated results file
    headers = [pivotParameter,variable]
    time_series = pd.read_csv(dirPath+filePrefix+"1.csv",names=headers)
    x1 = time_series[pivotParameter]
    x = x1[1::2] # Take every other value to match with hourly averages
    results = pd.DataFrame()
    results["Time"] = x

    # Get averaged original data
    od = pd.read_csv(dirPath+original_data)
    od = od[variable]
    od1 = [float(i) for i in od[0::2]]
    od2 = [float(j) for j in od[1::2]]
    od_avg = [(a+b)/2 for a,b in zip(od1,od2)]
    writetocsv = pd.DataFrame()
    writetocsv["Time"] = x
    writetocsv.insert(loc=1,column="Temperature",value=od_avg)
    writetocsv.to_csv(avg_orig_data,index=False)

    # Loop over all synthetic data output files to extract and write results
    for i in range(histories):
        df = pd.read_csv(dirPath+filePrefix+str(i)+".csv")
        y = df[variable]
        y1 = [float(j) for j in y[0::2]]
        y2 = [float(k) for k in y[1::2]]
        y_avg = [(a+b)/2 for a,b in zip(y1,y2)]
        results = results.rename(columns={'Data': 'Old_data'})
        results.insert(loc=1,column="Data",value=y_avg)

    results.to_csv(outputFile,index=False)

#elif plot_data == True:
    # Plot data on the same graph
data = pd.read_csv(outputFile,header=None, skiprows=1, index_col=0)
orig_data = pd.read_csv(avg_orig_data,header=None, skiprows=1, index_col=0)

#data = data.plot.line(title="ARMA Synthetic Temperature Histories",

#                        xlabel="Time (yrs)",
#                        ylabel=("Synthetic avergage hourly temperature (C)"),
#                        color="#929591",
#                        legend=None)
#orig_data = orig_data.plot.line(color="k",legend="Original temperature data")
#fig = data.get_figure()
#fig.savefig("temp_arma.png")
#fig = orig_data.get_figure()
#fig.savefig("temp_arma_2.png")





# Plot all columns of synthetic data in grey
ax = data.plot.line(title="ARMA Synthetic Temperature Histories",
                    xlabel="Time (yrs)",
                    ylabel=("Synthetic avergage hourly temperature (C)"),
                    color="#929591",
                    legend=None)

# Plot one column of original data in black
orig_data.rename(columns={orig_data.columns[0]: "Original temperature data"},
                 inplace=True)
orig_data.plot.line(color="k",
                    xlabel="Time (yrs)",
                    ax=ax)

ax.figure.show()
# Create and save figure
ax.figure.savefig("temp_arma.png")
