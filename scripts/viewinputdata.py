#! /usr/bin/env python

from runnerase import Runner
from runnerase import read, read_runnerconfig
from runnerase import generate_symmetryfunctions
from ase.visualize import view
import sys

if len(sys.argv)==1:
    runnerdata='input.data'
else:
    runnerdata=sys.argv[1]
dataset = read(runnerdata, index=':', format='runner')
view(dataset)
