#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler

# Load data from the file
data = np.loadtxt('energy.comp', skiprows=15)

Ha2eV = 27.211369917461
natoms = data[:, 1]
edft, enn = data[:, 2], data[:, 3]
diff = (edft - enn) / natoms

x = range(len(natoms))
ydft, ynn = edft * Ha2eV, enn * Ha2eV
ydiff = diff * Ha2eV

# Create plot
fig, ax1 = plt.figure(figsize=(6.75, 4.5), dpi=150), plt.subplot(1, 1, 1)
ax1.plot(x, ydft, color='tab:blue', label=r'$E_{DFT}$')
ax1.plot(x, ynn, color='tab:green', label=r'$E_{NN}$')

# Set labels and formatting
ax1.set_title('line plot with data points', fontsize=15)
ax1.set_xlabel(r'Configuration', fontsize=15)
ax1.set_ylabel(r'$E$ (eV)', fontsize=15)
ax1.tick_params(axis='both', labelsize=12)

# Create secondary y-axis
ax2 = ax1.twinx()
ax2.plot(x, ydiff, color='tab:orange', label='$E_{DFT}$ - $E_{NN}$')
ax2.set_ylabel(r'Difference (eV/atom)', fontsize=15)
ax2.spines['right'].set_color('tab:orange')
ax2.yaxis.label.set_color('tab:orange')
ax2.tick_params(axis='y', colors='tab:orange', labelsize=12)

# Combine legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0, prop={'size': 12})

plt.tight_layout()
plt.savefig("n2p2.png", dpi=300)
plt.show()
