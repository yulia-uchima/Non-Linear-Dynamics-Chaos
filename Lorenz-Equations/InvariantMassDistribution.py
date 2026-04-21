import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Values simulated by Monte Carlo-----------------------------------
# Load the simulated data for different processes related to the Higgs boson
monte_carlo = pd.read_csv("monte_carlo_higgs.csv")

# Convert relevant columns to NumPy arrays for easier plotting-----------------------------
zz = np.array(monte_carlo.zz)       # Array for ZZ events
dy = np.array(monte_carlo.dy)       # Array for DY events
ttbar = np.array(monte_carlo.ttbar) # Array for ttbar events
hzz = np.array(monte_carlo.hzz)     # Array for HZZ events

# Measured Data========================================================================

# Links to datasets containing measured data from CERN
links = [
    'http://opendata.cern.ch/record/5200/files/4mu_2011.csv',
    'http://opendata.cern.ch/record/5200/files/4e_2011.csv',
    'http://opendata.cern.ch/record/5200/files/2e2mu_2011.csv',
    'http://opendata.cern.ch/record/5200/files/4mu_2012.csv',
    'http://opendata.cern.ch/record/5200/files/4e_2012.csv',
    'http://opendata.cern.ch/record/5200/files/2e2mu_2012.csv'
]

# Load measured data from CSV files
df1 = pd.read_csv(links[0])
df2 = pd.read_csv(links[1])
df3 = pd.read_csv(links[2])
df4 = pd.read_csv(links[3])
df5 = pd.read_csv(links[4])
df6 = pd.read_csv(links[5])

# Concatenate all measured data into a single DataFrame
Data = pd.concat([df1, df2, df3, df4, df5, df6])

# Conditions for the region under consideration 
rmin = 70      # Minimum range for mass
rmax = 181     # Maximum range for mass
nbins = 37     # Number of bins for the histogram

# Create a histogram for the measured mass data within specified range
M_hist = np.histogram(Data['M'], bins=nbins, range=(rmin, rmax))

# Extract histogram values and bins
hist, bins = M_hist 
width = 1.0 * (bins[1] - bins[0])  # Calculate the width of each bin
center = (bins[:-1] + bins[1:]) / 2  # Calculate the center of each bin

# Histogram (Invariant mass distribution-)----------------------------------------------
plt.figure(figsize=(15, 10))  # Set the figure size

# ttbar---------------------------------------------------------------------------------
# Create a bar plot for ttbar events
ttbar_bar = plt.bar(center, ttbar, align='center', width=width, color='gray', 
                    linewidth=0, edgecolor='b', alpha=0.5, label=r'$t\bar{t}$')

# DY------------------------------------------------------------------------------------
# Create a bar plot for DY events, stacked on top of ttbar
dy_bar = plt.bar(center, dy, align='center', width=width, color='darkred', 
                  linewidth=0, edgecolor='black', alpha=0.5, bottom=ttbar, label='dy')

# ZZ----------------------------------------------------------------------------------------
# Create a bar plot for ZZ events, stacked on top of DY and ttbar
zz_bar = plt.bar(center, zz, align='center', width=width, color='teal', 
                  linewidth=0, edgecolor='black', alpha=0.5, bottom=ttbar+dy, label=r'ZZ')

# Measured data-------------------------------------------------------------------------------------------
# Plot the measured data points on top of the histogram
plt.plot(center, hist, linestyle='None', marker='o', color='black', markersize=5, label=r'Measured Data')

# Add titles and labels to the plot
plt.title('Histogram for the invariant mass reconstruction of the Higgs Boson', fontsize=12)
plt.xlabel('m [GeV]', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.ylim(0, 18)                # Set the y-axis limits
plt.xlim(rmin, rmax)          # Set the x-axis limits
plt.legend(fontsize=15)        # Add a legend to the plot
plt.savefig("./HiggsBoson.png")  # Save the figure as an image

# Hzz----------------------------
plt.figure(figsize=(10, 5))   # Create a new figure for the HZZ signal

xerrs = [width * 0.5 for i in range(0, nbins)]  # Calculate x-errors for error bars
yerrs = np.sqrt(hist)                             # Calculate y-errors based on histogram counts

# Create a bar plot for HZZ events, stacked
zz_bar = plt.bar(center, hzz, align='center', width=width, color='teal', 
                  linewidth=0, edgecolor='black', alpha=0.5, bottom=ttbar + dy, label=r'HZZ ')

# Measured data for HZZ plot
plt.plot(center, hist, linestyle='None', marker='o', color='black', markersize=5, label=r'Measured Data')

# Add titles and labels to the HZZ plot
plt.title('Considering the only Higgs signal Hzz', fontsize=15)
plt.xlabel('m [GeV]', fontsize=15)
plt.ylabel('Frequency', fontsize=15)
plt.ylim(0, 10)             # Set the y-axis limits for HZZ plot
plt.xlim(110, 140)         # Set the x-axis limits for HZZ plot
plt.legend(fontsize=15)     # Add a legend to the HZZ plot

plt.savefig("./hzz.png")    # Save the HZZ figure as an image