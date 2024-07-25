#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler
from sklearn.metrics import mean_squared_error
import statistics

# Load data from the file
data = np.loadtxt('diff.out', skiprows=15)

Ha2eV = 27.211369917461
x = data[:, 0]
natoms = data[:, 1]
edft = data[:, 4]
enn = data[:, 5]
fmean = data[:, 6]
fmax = data[:, 7]

diff = edft - enn

edftperatom = np.divide(edft, natoms)
ennperatom = np.divide(enn, natoms)
y = fmax * Ha2eV
y2 = fmean * Ha2eV
MEANmax = format(statistics.mean(y) * 1000, '.2f')
MEANmean = format(statistics.mean(y2) * 1000, '.2f')

# Create plot
fig, ax = plt.subplots(figsize=(5.85, 4.5), dpi=150)
ax.plot(x, y, 'o', markerfacecolor='none', color='tab:blue')
ax.plot(x, y2, '+', markerfacecolor='none', color='tab:red')

# Set labels and formatting
ax.set_title(r'$F_\mathrm{NN1} & F_\mathrm{NN2}$ difference', fontsize=15)
ax.set_xlabel(r'$configurations$', fontsize=15)
ax.set_ylabel(r'$F_\mathrm{NN}$ (eV/bohr)', fontsize=15)
ax.tick_params(axis='both', labelsize=12)
ax.ticklabel_format(useMathText=True)
plt.text(0.65, 0.85, f"Fmax MAE = {MEANmax} meV/bohr", horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=12, color='tab:blue')
plt.text(0.65, 0.15, f"F MAE = {MEANmean} meV/bohr", horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=12, color='r')
fig.tight_layout()

# Save and show plot
plt.savefig("n2p2-committee-FF.png", dpi=300)
plt.show()
