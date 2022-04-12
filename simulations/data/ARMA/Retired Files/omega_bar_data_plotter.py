import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft


df = pd.read_csv("C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r2/synData_2.csv")
#df = pd.read_csv("C:/Users/aidan/projects/neup-ies/simulations/data/ARMA/r2/hour_angle_parsed_data.csv")

#df.hist(cumulative=True, density=1)
#plt.show()

# Plot data
x = df["Normalised hour angle"]
y = df["Clearsky DNI"]
for i in range(1,len(y),1):
    if y[i] < 0:
        y[i] = 0
             
        
z2 = np.polyfit(x,y,6)
x1 = np.arange(-3500, 3500, step=12)
y1 = []
predict = np.poly1d(z2)
for n in x1:
    y1.append(predict(n))
plt.scatter(x,y)
plt.scatter(x1,y1)
plt.title("Clearsky DNI vs. second angle - synthetic data")
plt.xlabel("Normalised second angle")
plt.ylabel("Clearsky DNI")
#plt.xticks([])
plt.xticks(np.arange(-3600, 3601, step=1200))
#plt.xticks(np.arange(-1.1, 1, step=0.2))
plt.show()

z2 = np.polyfit(x,y,6)

print(z2)

"""fy = fft(y1)
N = len(fy)

n = np.arange(N)
# get the sampling rate
sr = 1 / (12)
T = N/sr
freq = n/T 

# Get the one-sided specturm
n_oneside = N//2
# get the one side frequency
f_oneside = freq[:n_oneside]

plt.figure(figsize = (12, 6))
plt.plot(f_oneside, np.abs(fy[:n_oneside]), 'b')
plt.xlim(0, 1/1400)
t_h = 1/f_oneside / (1200)
plt.figure(figsize=(12,6))
plt.plot(t_h, np.abs(fy[:n_oneside])/n_oneside)
plt.xticks([1,1.27,1.82,2,2.18,6,4])
plt.xlim(0, 10)
plt.xlabel('Period ($hour angles$)')
plt.show()"""