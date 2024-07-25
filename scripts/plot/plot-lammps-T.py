#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler
import matplotlib.ticker as mtick

# Load data from the file
data = np.loadtxt('data.txt', skiprows=5, usecols=(2, 3))
x, y = data[:, 0], data[:, 1]

# Calculate mean value for y after 0.2 ns
mean_value = np.mean(y[66:])
mean_str = format(mean_value, '.2f')

# Create plot
fig, ax1 = plt.figure(figsize=(6.75, 4.5), dpi=150), plt.subplot(1, 1, 1)
ax1.plot(x, y, color='tab:blue', label=r'$Density_{mean}$')

# Set labels and formatting
ax1.set_xlabel(r'time/ps', fontsize=15)
ax1.set_ylabel(r'Temperature/K', fontsize=15)
ax1.tick_params(axis='both', labelsize=12)
plt.text(0.75, 0.15, f"MEAN = {mean_str} (after 0.2 ns)", horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes, fontsize=12, color='r')
plt.tight_layout()

# Save and show plot
plt.savefig("T.png", dpi=300)
plt.show()
print(mean_str)
