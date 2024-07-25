#! /usr/bin/env python
from ase.io import read

# Read POSCAR file
pos = read('POSCAR')
frompos = pos.cell.cellpar()

# Ensure the cell dimensions do not exceed 40 in any direction
if frompos[2] > 40:
    pos.set_cell([frompos[0], frompos[1], 40, frompos[3], frompos[4], frompos[5]])
    print("Adjusted z dimension to 40")
if frompos[1] > 40:
    pos.set_cell([frompos[0], 40, frompos[2], frompos[3], frompos[4], frompos[5]])
    print("Adjusted y dimension to 40")
if frompos[0] > 40:
    pos.set_cell([40, frompos[1], frompos[2], frompos[3], frompos[4], frompos[5]])
    print("Adjusted x dimension to 40")

# Center and wrap the positions
pos.wrap()
pos.center()
pos.set_center_of_mass((0.5, 0.5, 0.5), scaled=True)
pos.write('POSCAR')

# Check if the slab is breaking
if frompos[2] == 40:
    pos.set_cell([frompos[0], frompos[1], 20, 90, 90, 90])
    pos.wrap()
    pos.center()
    pos.set_center_of_mass((0.5, 0.5, 0.5), scaled=True)
    z_positions = pos.get_positions()[:, 2]
    if any(z < 0 for z in z_positions):
        pos.write("POSCARerrorinz")
        print("POSCAR error in z")
