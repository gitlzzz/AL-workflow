#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler
from matplotlib import font_manager
import os

# Load custom font
font_path = os.path.expanduser('~/soft/anaconda3/fonts/arial.ttf')
font_prop = font_manager.FontProperties(fname=font_path)

# Load data from the file
data = np.loadtxt('data.txt', skiprows=3, usecols=(2, 4))
x, y = data[:, 0], data[:, 1]

# Calculate mean value for y over the last 267 points
mean_value = np.mean(y[-267:])
mean_str = format(mean_value, '.4f')

# Create plot
fig, ax1 = plt.figure(figsize=(6.75, 4.5), dpi=150), plt.subplot(1, 1, 1)
ax1.plot(x, y, color='tab:blue', label=r'$Density_{mean}$')

# Set y-axis limits with margin
y_margin = 0
y_min, y_max = np.min(y) - y_margin, np.max(y) + y_margin
ax1.set_ylim(y_min, y_max)

# Set labels and formatting
ax1.set_xlabel(r'time/ps', fontproperties=font_prop, fontsize=15)
ax1.set_ylabel(r'potential energy/eV', fontproperties=font_prop, fontsize=15)
ax1.tick_params(axis='both', labelsize=12)
plt.text(0.75, 0.15, f"MEAN = {mean_str} (last 0.8 ns)", horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes, fontproperties=font_prop, fontsize=12, color='r')
plt.tight_layout()

# Save and show plot
plt.savefig("E.png", dpi=300)
# plt.show()
print(mean_str)
