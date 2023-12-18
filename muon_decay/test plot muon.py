import os
import fnmatch
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import numpy as np
import sys
sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\utils')
from parameters import *
current_directory = os.getcwd()
data_directory = current_directory+"\\muon_decay\\extracted_data\\"

files = ["Run_13_full_analysed\\time_intervals_run_000013.txt"]

if len(files) == 0:
    files = fnmatch.filter(os.listdir(data_directory), 'time_interval*')
time_intervals = np.array([])
for file in files:
    ti = np.loadtxt(data_directory+file)
    time_intervals = np.append(time_intervals,ti)
time_intervals*=clockcycle_value
time_intervals.sort()

# Exponential function for the fit
def exponential(x,A,tau):
    return A*np.exp(-x/tau)

from scipy.optimize import curve_fit

fig, ax = plt.subplots(1, 1, figsize=(6, 4))

# Create histogram
n, bins, _ = ax.hist(time_intervals, 100, density=False, alpha=0.7, label='Histogram')

# Get bin centers and bin width
bin_centers = (bins[1:] + bins[:-1]) / 2
bin_width = bins[1] - bins[0]

# Fit the exponential function to the data
coeffs, cov = curve_fit(exponential, bin_centers, n, p0=[1700, 2000])

# Plot the exponential fit
t = np.linspace(min(time_intervals), max(time_intervals), 100)
ax.plot(t, exponential(t, *coeffs), '-r', lw=3, label='Exponential Fit')

# Add error bars to the histogram bars
errors = np.sqrt(n)  # Simple example, you may want to customize based on your data
ax.errorbar(bin_centers, n, yerr=errors, fmt='none', color='black', capsize=3, label='Error Bars')

# Set labels and legend
ax.set_xlabel("$\Delta t$")
ax.set_ylabel("$n_{events}$")
ax.set_ylim(0, max(n)*1.1)
ax.legend()

plt.show()
fig.savefig("muon_lifetime.pdf", format="pdf", bbox_inches="tight")

tau = coeffs[1]
tau_std = np.sqrt(cov[1,1])


print("Muon lifetime : \u03C4 = {} \u00B1 {} ns".format(tau,tau_std))