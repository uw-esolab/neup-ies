import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import numpy as np
import statistics as stats

filename_orig = "r1/solar_data_parsed.csv"

headers = ["Time","DHI"]

original_data = pd.read_csv(filename_orig,names=headers)

x = original_data["Time"]
y = original_data["DHI"]

xa = np.array(x[1:])
xb = xa.astype(float)
ya = np.array(y[1:])
yb = ya.astype(float)


plt.plot(xb,yb)
plt.xticks(np.arange(0, len(xb)*30, 100000))
plt.yticks(np.arange(0, 500, 100))
plt.show()

fy = fft(yb)
N = len(fy)

n = np.arange(N)
# get the sampling rate
sr = 1 / (30*60)
T = N/sr
freq = n/T 

# Get the one-sided specturm
n_oneside = N//800
# get the one side frequency
f_oneside = freq[:n_oneside]

plt.figure(figsize = (12, 6))
plt.plot(f_oneside, np.abs(fy[:n_oneside]), 'b')
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|')
plt.show()


plt.plot(xb,fy)
plt.title("Fast Fourier transform of original solar data")
plt.xlabel("Frequency (hz)")
plt.ylabel("Amplitude")
#plt.xticks(np.arange(0, len(xb), 10))
plt.show()
