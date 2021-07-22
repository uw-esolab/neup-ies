import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import numpy as np
import statistics as stats

filename_orig = "r1/solar_data_parsed.csv"

headers = ["Time","DNI"]

original_data = pd.read_csv(filename_orig,names=headers)

x = original_data["Time"]
y = original_data["DNI"]

y1 = [float(i) for i in y[1::2]]
y2 = [float(j) for j in y[2::2]]

y_avg = [(a+b)/2 for a,b in zip(y1,y2)]

x2 = x[1::2]

y = [float(i) for i in y[1:]]
x = x[1:]
a = fft(y)
plt.plot(x,a)
plt.title("Fast Fourier transform of original solar data")
plt.xlabel("Frequency (hrs)")
plt.ylabel("Amplitude")
plt.xticks(np.arange(0, len(x), 1000))
plt.show()
