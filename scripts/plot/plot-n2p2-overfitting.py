#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler
import matplotlib.ticker as mtick

# Load data from the file
data = np.loadtxt('learning-curve.out', skiprows=34)

Ha2eV = 27.211386245988  # RuNNerUC
x = data[:, 0]
Etrain = data[:, 1] * Ha2eV
Etest = data[:, 2] * Ha2eV
Ftrain = data[:, 9] * Ha2eV
Ftest = data[:, 10] * Ha2eV

# Create plot
fig, ax1 = plt.subplots(figsize=(6.75, 4.5), dpi=150)
ax1.plot(x, Etrain, color='tab:blue', label=r'$E_{train}$')
ax1.plot(x, Etest, color='tab:blue', linestyle='dashed', label=r'$E_{test}$')
ax1.set_yscale('log')
ax1.set_xlabel(r'epoch', fontsize=15)
ax1.set_ylabel(r'$E$ (eV/atom)', fontsize=15)
ax1.tick_params(axis='both', labelsize=12)
ax1.yaxis.label.set_color('tab:blue')
ax1.spines['right'].set_color('tab:blue')
ax1.tick_params(axis='y', colors='tab:blue', labelsize=12)

# Create secondary y-axis
ax2 = ax1.twinx()
ax2.plot(x, Ftrain, color='tab:orange', label='$F_{train}$')
ax2.plot(x, Ftest, color='tab:orange', linestyle='dashed', label='$F_{test}$')
ax2.set_yscale('log')
ax2.set_ylabel(r'Force (eV/bohr)', fontsize=15)
ax2.spines['right'].set_color('tab:orange')
ax2.yaxis.label.set_color('tab:orange')
ax2.tick_params(axis='y', colors='tab:orange', labelsize=12)

# Combine legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0, prop={'size': 12})

plt.tight_layout()
plt.savefig("n2p2-fiting.png", dpi=300)
plt.show()
