#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler

# Load data from the file
data = np.loadtxt('diff.out', skiprows=1)

Ha2eV = 27.211369917461
x = data[:, 0]
natoms = data[:, 1]
enn1, enn2 = data[:, 4], data[:, 5]
diff = enn1 - enn2
ynn1, ynn2 = enn1 * Ha2eV, enn2 * Ha2eV
ydiff = diff * Ha2eV / natoms

# Create plot
fig, ax1 = plt.figure(figsize=(9, 6), dpi=150), plt.subplot(1, 1, 1)
ax1.plot(x, ynn1, color='tab:blue', label=r'$E_{NN1}$')
ax1.plot(x, ynn2, color='tab:green', label=r'$E_{NN2}$')

# Set labels and formatting
ax1.set_title('line plot with data points', fontsize=20)
ax1.set_xlabel(r'Configuration', fontsize=20)
ax1.set_ylabel(r'$E$ (eV)', fontsize=20)
ax1.tick_params(axis='both', labelsize=15)

# Create secondary y-axis
ax2 = ax1.twinx()
ax2.plot(x, ydiff, color='tab:orange', label='$E_{NN1}$ - $E_{NN2}$')
ax2.set_ylabel(r'Difference (eV/atom)', fontsize=20)
ax2.spines['right'].set_color('tab:orange')
ax2.yaxis.label.set_color('tab:orange')
ax2.tick_params(axis='y', colors='tab:orange', labelsize=15)

# Combine legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0, prop={'size': 15})

plt.tight_layout()
plt.savefig("n2p2-committee.png", dpi=300)
plt.show()
