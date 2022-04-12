import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import numpy as np
import statistics as stats

filename_orig = "r2/hour_angle_parsed_data_edit.csv"

#headers = ["Normalised hour angle","Clearsky DNI"]

original_data = pd.read_csv(filename_orig)

x = original_data["Normalised hour angle"]
y = original_data["Clearsky DNI"]

xa = np.array(x[1:])
xb = xa.astype(float)
ya = np.array(y[1:])
yb = ya.astype(float)


plt.scatter(xb,yb)
plt.title("Clearsky DNI vs. hour angle - synthetic data")
plt.xlabel("Normalised hour angle")
plt.ylabel("Clearsky DNI")
#plt.xticks([])
plt.xticks(np.arange(-60, 60, step=10))
plt.show()

fy = fft(yb)
N = len(fy)

n = np.arange(N)
# get the sampling rate
sr = 1 / (30*60)
T = N/sr
freq = n/T 

# Get the one-sided specturm
n_oneside = N//2
# get the one side frequency
f_oneside = freq[:n_oneside]

plt.figure(figsize = (12, 6))
plt.plot(f_oneside, np.abs(fy[:n_oneside]), 'b')
t_h = 1/f_oneside / (3600)
plt.figure(figsize=(12,6))
plt.plot(t_h, np.abs(fy[:n_oneside])/n_oneside)
plt.xticks([1,1.27,1.82,2,2.18,6,4])
plt.xlim(0, 10)
plt.xlabel('Period ($hour angles$)')
plt.show()
