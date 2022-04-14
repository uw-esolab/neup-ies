import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq
import numpy as np
import statistics as stats
from peakdetect import peakdetect


filename_orig = "r1/synData_8.csv"
#filename_orig = "r1/solar_data_without_clearsky.csv"
filename_orig_time = "r2/hour_angle_parsed_data_reverse.csv"
filename_orig_DNI = "r2/synData_2.csv"


#headers1 = ["Time1"	,"Time2","Clearsky DNI",	"Hour angle","Sunrise hour angle",	"Sunset hour angle",	"Normalised hour angle","clear1","clear2", "New Normalised hour second",	"New Clearsky DNI",	"New Normailised hour angle"]
headers = ["Time", "DNI"]
#headers = ["Time1"	,"Time","Clearsky DNI",	"Hour angle","Sunrise hour angle",	"Sunset hour angle",	"Normalised hour angle"]
synheaders = ["Normalised hour angle", "Clearsky_DNI"]

original_data_time = pd.read_csv(filename_orig,names=headers)
original_data_DNI = pd.read_csv(filename_orig ,names=headers)

x = original_data_time["Time"]
y = original_data_DNI["DNI"]

xa = np.array(x[1:])
xb = xa.astype(float)
ya = np.array(y[1:])
yb = ya.astype(float)

"""for i in range(len(yb)):
    if yb[i] < 600:
        yb[i] = 0"""

plt.plot(xb,yb)
#plt.xlim(0, 25000)
plt.xticks(np.arange(0, 600000, 100000))
plt.yticks(np.arange(0, 1100, 100))
plt.show()

plt.plot(xb,yb)
plt.xlim(0, 25000)
#plt.xticks(np.arange(0, 600000, 100000))
#plt.yticks(np.arange(0, 1100, 100))
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
plt.xlabel('Frequency ($Hz$)')
plt.loglog(f_oneside, np.abs(fy[:n_oneside]), 'b')
t_h = 1/f_oneside / (60 * 60)
#print(type(t_h))
peaks  = peakdetect(t_h)
print(peaks)
plt.figure(figsize=(12,6))
plt.plot(t_h, np.abs(fy[:n_oneside])/n_oneside)
plt.xticks([24,9.6,12,8,6,4.8,4,(24/7),3, 24/9, 24/10, 24/10, 24/12])
plt.xlim(0, 50)
plt.xlabel('Period ($hour$)')
plt.show()


"""# FFT the signal
sig_fft = fft(yb)
# copy the FFT results
sig_fft_filtered = sig_fft.copy()

N = len(fy)

n = np.arange(N)
# get the sampling rate
sr = 1 / (60*60)
T = N/sr
freq = n/T 

# define the cut-off frequency
cut_off = 0.000125

# high-pass filter by assign zeros to the 
# FFT amplitudes where the absolute 
# frequencies smaller than the cut-off 
sig_fft_filtered[np.abs(freq) < cut_off] = 0

# get the filtered signal in time domain
filtered = ifft(sig_fft_filtered)

# plot the filtered signal
plt.figure(figsize = (12, 6))
plt.plot(xb, filtered)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.show()

# plot the FFT amplitude before and after
plt.figure(figsize = (12, 6))
plt.subplot(121)
plt.stem(freq, np.abs(sig_fft), 'b', \
         markerfmt=" ", basefmt="-b")
plt.title('Before filtering')
plt.xlim(0, 0.00005)
plt.xlabel('Frequency (Hz)')
plt.ylabel('FFT Amplitude')
plt.subplot(122)
plt.stem(freq, np.abs(sig_fft_filtered), 'b', \
         markerfmt=" ", basefmt="-b")
plt.title('After filtering')
plt.xlim(0.0005, 0.0006)
plt.xlabel('Frequency (Hz)')
plt.ylabel('FFT Amplitude')
plt.xticks([0.00054375])
plt.tight_layout()
plt.show()"""