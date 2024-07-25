#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler
from sklearn.metrics import mean_squared_error

# Load data from the file
data = np.loadtxt('energy.comp', skiprows=15)

Ha2eV = 27.211369917461
natoms = data[:, 1]
edft, enn = data[:, 2], data[:, 3]
diff = edft - enn

edftperatom = np.divide(edft, natoms)
ennperatom = np.divide(enn, natoms)
x, y = edftperatom * Ha2eV, ennperatom * Ha2eV
RMSE = format(mean_squared_error(x, y, squared=False) * 1000, '.2f')

# Create plot
fig, ax = plt.subplots(figsize=(5.85, 4.5), dpi=150)
ax.plot(x, x, color='black')
ax.plot(x, y, 'o', markerfacecolor='none', color='tab:blue')

# Set labels and formatting
ax.set_title(r'$E_\mathrm{DFT} & E_\mathrm{NN}$', fontsize=15)
ax.set_xlabel(r'$E_\mathrm{DFT}$ (eV/atom)', fontsize=15)
ax.set_ylabel(r'$E_\mathrm{NN}$ (eV/atom)', fontsize=15)
ax.tick_params(axis='both', labelsize=12)
ax.ticklabel_format(useMathText=True)
plt.text(0.75, 0.15, f"RMSE = {RMSE} meV/atom", horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes, fontsize=12, color='r')
fig.tight_layout()

# Save and show plot
plt.savefig("n2p2-ee.png", dpi=300)
plt.show()
