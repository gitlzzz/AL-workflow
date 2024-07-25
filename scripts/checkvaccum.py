#! /usr/bin/env python
###fail to detect vaccum direction when structure is highly tri
import numpy as np
from ase.io import read, write
from pathlib import Path

def bottom(pos, axis):
    coord=pos.get_positions()[:,axis]
    bottom = min(coord)
    pos.positions[:,axis]=coord - bottom + 2
    pos.write('POSCAR')

pos=read('POSCAR')
pos.wrap() # to make sure atoms in the same cell
#pos=pos*(2,2,2)
pos.wrap() 
#cellsize=pos.cell.cellpar()
cell = pos.cell
cellsize=[cell[0][0],cell[1][1],cell[2][2]]
print(cellsize)
allz=pos.get_positions()

#sort the value of coordinates along each axis
Xcoord=np.array(allz[:,0])
Xcoord=np.sort(Xcoord)
print(Xcoord)
Ycoord=np.array(allz[:,1])
Ycoord=np.sort(Ycoord)
Zcoord=np.array(allz[:,2])
Zcoord=np.sort(Zcoord)

#calculate the maximum difference as vaccum layer thickness
num=len(Xcoord)
Xdiff=np.array([None for x in range(num)])
Ydiff=np.array([None for x in range(num)])
Zdiff=np.array([None for x in range(num)])
for i in range(0,num-1):
    Xdiff[i]=Xcoord[i+1]-Xcoord[i]
    Ydiff[i]=Ycoord[i+1]-Ycoord[i]
    Zdiff[i]=Zcoord[i+1]-Zcoord[i]
Xdiff[-1]=Xcoord[0]+cellsize[0]-Xcoord[-1]
Ydiff[-1]=Ycoord[0]+cellsize[1]-Ycoord[-1]
Zdiff[-1]=Zcoord[0]+cellsize[2]-Zcoord[-1]

alldiff=np.concatenate((Xdiff,Ydiff,Zdiff))
print(Xdiff)
print(Ydiff)
print(Zdiff)
#print(alldiff)
maxalongaxis=np.array([None for x in range(3)])
maxalongaxis[0]=np.amax(Xdiff)
maxalongaxis[1]=np.amax(Ydiff)
maxalongaxis[2]=np.amax(Zdiff)

#check if the vaccum layer is less for VDW and just a little bit more than bonding.
if ((10.0 > alldiff) & (alldiff > 3.0)).any():
    pos.write('POSCARerrorinz')
    if ((10.0 > alldiff) & (alldiff > 6.0)).any():
        pos.write('POSCARneedtune')
elif (maxalongaxis >= 10.0).all():
    Path('vaccumalongxyz').touch()
    bottom(pos,2)
elif maxalongaxis[0] >= 10.0:
    Path('vaccumalongx').touch()
    bottom(pos,0)
elif maxalongaxis[1] >= 10.0:
    Path('vaccumalongy').touch()
    bottom(pos,1)
elif maxalongaxis[2] >= 10.0:
    Path('vaccumalongz').touch()
    bottom(pos,2)
