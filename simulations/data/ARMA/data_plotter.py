import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import numpy as np
import statistics as stats

# Read all files
filename_orig = "r1/temperature_data_parsed.csv"
filename_arma = "r1/temperature_data_arma.csv"
filename_farma = "r1/temperature_data_farma.csv"

headers = ["Time","Temperature"]

original_data = pd.read_csv(filename_orig,names=headers)
arma_data = pd.read_csv(filename_arma,names=headers)
farma_data =pd.read_csv(filename_farma,names=headers)

# Extract two columns of data
x = original_data["Time"]
y = original_data["Temperature"]
z = arma_data["Temperature"]
w = farma_data["Temperature"]

# Get averages per hour
y1 = [float(i) for i in y[1::2]]
y2 = [float(j) for j in y[2::2]]

z1 = [float(i) for i in z[1::2]]
z2 = [float(j) for j in z[2::2]]

w1 = [float(i) for i in w[1::2]]
w2 = [float(j) for j in w[2::2]]

y_avg = [(a+b)/2 for a,b in zip(y1,y2)]
z_avg = [(a+b)/2 for a,b in zip(z1,z2)]
w_avg = [(a+b)/2 for a,b in zip(w1,w2)]

# Calculate statistical parameters
y_stdev = stats.stdev(y_avg)
y_text = "Stdev = " + str("{:.2f}".format(y_stdev))
z_stdev = stats.stdev(z_avg)
z_text = "Stdev = " + str("{:.2f}".format(z_stdev))
w_stdev = stats.stdev(w_avg)
w_text = "Stdev = " + str("{:.2f}".format(w_stdev))

x2 = x[1::2]

a = fft(y1)
plt.plot(x2,a)
plt.title("Fast Fourier transform of original temperature data")
plt.xlabel("Frequency (hrs)")
plt.ylabel("Amplitude")
plt.xticks(np.arange(0, len(x2), 1000))
plt.show()

plt.subplot(3,1,1)
plt.plot(x2,y_avg)
plt.title("Original input temperature data")
plt.yticks(np.arange(-30, 50, 10))
plt.xticks(np.arange(0, len(x2), 1000))
plt.xlabel("Time (hrs)")
plt.ylabel("Average hourly temperature (C)")
plt.text(8700,40, y_text)

plt.subplot(3,1,2)
plt.plot(x2,z_avg)
plt.title("ARMA output temperature data")
plt.yticks(np.arange(-30, 50, 10))
plt.xticks(np.arange(0, len(x2), 1000))
plt.xlabel("Time (hrs)")
plt.ylabel("Synthetic avg hourly temperature (C)")
plt.text(8700,40, z_text)

plt.subplot(3,1,3)
plt.plot(x2,w_avg)
plt.title("FARMA output temperature data")
plt.yticks(np.arange(-30, 50, 10))
plt.xticks(np.arange(0, len(x2), 1000))
plt.xlabel("Time (hrs)")
plt.ylabel("Synthetic avg hourly temperature (C)")
plt.text(8700,40, w_text)

plt.subplots_adjust(hspace=0.5)
plt.show()
