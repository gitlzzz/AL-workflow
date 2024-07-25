#!/usr/bin/python

'''
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\                                                      /\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/                 RuNNerActiveLearn.py                 \/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\                  Version 02-07-2021                  /\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/                    Marco Eckhoff                     \/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\         Georg-August-Universitaet Goettingen         /\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/   Theoretische Chemie, Institut fuer Physikalische   \/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\ Chemie, Tammannstrasse 6, 37077 Goettingen, Germany  /\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/        marco.eckhoff@chemie.uni-goettingen.de        \/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\                                                      /\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
/\                                                                                                /\ 
\/ RuNNerActiveLearn sets up and analyses LAMMPS simulations to find structures exhibiting inter- \/ 
/\ and extrapolation errors of the underlying high-dimensional neural network potential (HDNNP).  /\ 
\/ The aim of the algorithm is to create a selection of uncorrelated structures which are missing \/ 
/\ in the current reference data set of the HDNNP. These structures are written to the            /\ 
\/ input.data-add file and have to be recalculated by the reference method to improve the         \/ 
/\ reliability of the HDNNP. RuNNerActiveLearn automatically adjusts some required thresholds in  /\ 
\/ this self-learning procedure and provides reasonable recommendations for the others to sample  \/ 
/\ the configuration space efficiently. Interpolation errors are detected by a comparison of the  /\ 
\/ results of two different HDNNPs. Statistics of the extrapolated structures are provided as     \/ 
/\ well. The program package requires binaries of LAMMPS including the n2p2 libraries and RuNNer. /\ 
\/                                                                                                \/ 
S/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
/\                                                                                                /\ 
\/ For using RuNNerActiveLearn please cite:                                                       \/ 
/\ M. Eckhoff and J. Behler, arXiv:2104.14439 [physics.comp-ph], (2021).                          /\ 
\/ M. Eckhoff and J. Behler, J. Chem. Theory Comput. 15, 3793-3809 (2019).                        \/ 
/\                                                                                                /\ 
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
\/                                                                                                \/ 
/\ RuNNerActiveLearn has been developed and written by Marco Eckhoff.                             /\ 
\/                                                                                                \/ 
/\ Contributions were made by                                                                     /\ 
\/ Joerg Behler: General procedure                                                                \/ 
/\ Alea Miako Tokita: Testing                                                                     /\ 
\/ Knut Nikolas Lausch: Testing                                                                   \/ 
/\                                                                                                /\ 
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
\/                                                                                                \/ 
/\ This program is free software: you can redistribute it and/or modify it under the terms of the /\ 
\/ GNU General Public License as published by the Free Software Foundation, either version 3 of   \/ 
/\ the License, or (at your option) any later version.                                            /\ 
\/ This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;      \/ 
/\ without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  /\ 
\/ the GNU General Public License for more details.                                               \/ 
/\ You should have received a copy of the GNU General Public License along with this program. If  /\ 
\/ not, see http://www.gnu.org/licenses.                                                          \/ 
/\                                                                                                /\ 
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ 
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
'''

####################################################################################################

####################################################################################################

import numpy as np
import sys
import os
import subprocess
import warnings
from copy import deepcopy

####################################################################################################

####################################################################################################

'''
Program input: Adjust the following settings for your system.

'''

# Set an integer to specify the mode. Mode 1 generates LAMMPS input files and mode 2 compiles the
# results and selects the important new reference structures. The setting can be overwritten by a
# command line argument.
mode = 1

# Set an integer to define the seed = 1
# 1. If you run RuNNerActiveLearn for another cycle of including missing structures, set this value
# to a higher value than the highest value of the last cycle in order to be sure that every LAMMPS
# simulation is different.
seed = 1

# Set an array of strings which defines the usage of the nve, nvt, and/or npt integrators.
# Recommendation: Use the npt integrator as this varies the density/simulation cell as well.
integrators = ['npt']
# Set an array of integers/floats of temperature values in K.
# Recommendation: Set the upper limit to a value 25 to 50% higher than the desired highest
# temperature in later simulations if the system does not undergo any phase transition at this
# temperature to ensure that the HDNNP is well sampled at the desired highest temperature.
temperatures = range(300, 801, 50)
# Set an array of integers/floats of pressure values in bar used in NpT simulations. Set an empty
# array if no NpT simulations are performed.
# Recommendation: Use [1] for solids as this normally does not have to be varied.
pressures = [1]
# Set an integer to define the number of MD steps (simulation.lammps has to use the variable N_steps
# as well).
# Recommendation: N_steps * timestep = 100 ps.
N_steps = 100000 
# Set a float to define the timestep in ps.
# Recommendation: 0.0005 for systems with H atoms, otherwise 0.001.
timestep = 0.003
# Set a string which specifies the barostat option of the npt integrator: iso, aniso, or tri (iso
# and aniso are not supported in combination with a non orthorhombic cell). If the npt integrator is
# not used, choose the option according to the simulation cell. Non orthorhombic cells require to
# set tri otherwise iso or aniso can be selected (no difference if npt integrator is not used).
# Recommendation: Set to tri as this is without any restrictions. But dependent on the system there
# can be good reasons for the usage of one of the others.
barostat_option = 'tri'
# Set a string which specifies the atom style of LAMMPS structure.lammps file: atomic or full.
# Recommendation: full is currently only required for magnetic HDNNPs (requires to compile LAMMPS
# with the molecule package).
atom_style = 'atomic'
# Set an integer which defines that only every nth structure is kept in the structures.lammpstrj
# file if no extrapolation occured. The value has to be a divisor of N_steps.
# Recommendation: dump_lammpstrj * timestep = 0.1 ps.
dump_lammpstrj = 40
# Set an integer/float to define the RuNNer cutoff radius in Bohr radii.
# Recommendation: 12.0.
runner_cutoff = 12.0
# Set to True for periodic systems and to False for non-periodic systems. For non-periodic systems
# an orthorhombic simulation cell with lattice constants (x_max - x_min + 2 * runner_cutoff, y_max
# - y_min + 2 * runner_cutoff, z_max - z_min + 2 * runner_cutoff) is used in the LAMMPS simulation.
periodic = True
# Set an array of strings which defines the element symbols in elemental order.
element_types = ['O', 'Cu']
# Set an array of integer/floats which defines the element masses in atomic units in elemental
# order.
masses = [15.99491462, 63.546]
# Set an integer which specifies the maximal length of the job list including the LAMMPS simulations
# as there might be a limit of the job array size. If the number of jobs is higher, several job
# lists are generated. If the value is set to 0, all jobs are compiled in one job list.
max_len_joblist = 280 #356 on marenostrum

# Set a string which specifies the line start of lines in the input.data file which include the name
# of the structure. This enables to apply different settings for different groups of structures, for
# example, if their training progress is at different levels. Furthermore, the name might be
# required for the assignment of settings in the later electronic structure calculations. Avoid '_'
# in this structure name as this sign is used in the following as separator. If a structure name is
# not required, None has to be assigned. If None is assigned, comment_name_index and
# comment_name_separator will not be used.
#one structure one 'comment'
comment_name_keyword = 'comment'
# Set an integer which specifies the index of the name in the above given lines (counting starts at
# 0).
comment_name_index = 1
# Set a string which specifies an additional separator of the name. Only the substring before this
# separator is used.
# Recommendation: '_'
comment_name_separator = ' '
# Set an array of strings of used structure names or set None if they shall not be used.
structure_names = ['cuxobulk','cuxosurf']
# If different structure names are used, the following settings can be specified for each structure
# name individually or one value has to be specified which is used for all structure names. In any
# case the following settings require the specification of an array. The order has to be the same as
# for the structure names. If no structure names are used or if a value shall be used for all
# structures, arrays with one entry are required.

# Set an array of arrays of selections to use only every jth structure starting with the ith one for
# each structure name (format: [[i, j], ...]).
# Recommendation: This depends on the reliability of the HDNNP. The usage the current input.data
# file as source of new structures is recommented. You should not add more than about a third of the
# current number of structures, i.e., performing too many simulations is not efficient. The more
# often you retrain your HDNNP, the less likely it is too include structures fixing the same
# problem. However, if you include too less structures in every iteration, the procedure will also
# be slow. If you try to find the last gaps in the sampled configuration space, you can do more
# simulations and also longer simulations.
structure_selection = [[int(seed),20],[int(seed),20]]

# Set an array of dictionaries of minimal nearest neighbour distances for each structure name
# (format: [{A:[A-A, A-B, A-C], B:[B-B, B-C], C:[C-C]}, ...] where A, B, and C are strings of the
# element symbols in the elemental order A, B, C and A-A, A-B, ..., and C-C are floats of the
# minimal interatomic distances in Angstrom).
d_mins = [{'O':[1.07, 1.51], 'Cu':[1.83]}]
# Set an array of integers to define the minimal number of time steps between two extrapolated
# structures for each structure name.
# Recommendation: min_timestep_separation_extrapolation * timestep = 0.01 ps.
min_timestep_separation_extrapolation = [3,3]
# Set an array of integers to define the usual time step separations between two structures for
# interpolation checks for each structure name (this can be smaller in case only less than three
# checks would be possible). The value has to be smaller than a fifth of N_steps.
# Recommendation: timestep_separation_interpolation_checks * timestep = 5 ps.
timestep_separation_interpolation_checks = [1000,1000]
# Set an array of integers to define the minimal number of time steps between two structures for
# interpolation checks for each structure name (this can lead to less than three checks).
# Recommendation: min_timestep_separation_interpolation * timestep = 0.1 ps.
min_timestep_separation_interpolation = [40,40]
# Set an array of floats to define the thresholds for the maximal energy difference between both
# HDNNP predictions in Hartree/atom for each structure name.
# Recommendation: Highest energy error of the training data set (excluding outlyers) or five times
# the energy RMSE of the training dat set for the given structures.
delta_E = [0.00090,0.00090]
# Set an array of floats to define the threshold for the maximal force component difference between
# both HDNNP predictions in Hartree/Bohr radius for each structure name.
# Recommendation: 100 times the value of delta_E or highest force error of the training data set
# (excluding outlyers) or five times the force RMSE of the training dat set for the given
# structures.
delta_F = [0.00640,0.00640]
# Set an array of Boolean values which specifies for each structure name if all extrapolated
# structures shall be selected. Otherwise they are only selected if they above the energy and force
# thresholds. max_extrapolated_structures and exceptions specified below can deselect the
# extrapolated structures again.
# Recommendation: Set to True if the HDNNP has reached some reliability.
all_extrapolated_structures = [True,True]
# Set an array of integers which defines the maximal number of the same kind of extrapolations for
# each structure name. Set to 0 if no limit should be applied.
# Recommendation: In the initial generations of improvements it can happen that a particular
# symmetry function leads to most of the extrapolations. In this case it makes sense to select only
# a fraction of these structures since they fix all the same problem. A value about 50 might be
# reasonable.
max_extrapolated_structures = [500,500]
# Set an array of integers which defines the usual maximal numbers of selected interpolated
# structures per simulation for each structure name.
# Recommendation: [4].
max_interpolated_structures_per_simulation = [50,50]
# Set an array of arrays/None to manually define exceptions of small extrapolations which shall be
# only included to a certain fraction for each structure name. The max_extrapolated_structures
# limitation is overwritten for the given extrapolations. An example for the format is
# [[[A, B, C], ...], None, ...] where A is a string of the element symbols specifying the central
# atoms of the extrapolated symmetry functions, B is a string of the numbers of corresponding
# symmetry functions, and C is a float of the used fraction. A and B have to be identical to the
# entries in input.data-new. For each structure name an array of several exceptions ([A, B, C]) can

# input.data-new and extrapolation_statistics_XX.dat.
exceptions = [None,None]

# Set a string which is used as command to run the RuNNer executable.
#RuNNer_exe = 'RuNNer.x'
RuNNer_exe = 'RuNNer-3-gfortran.x'
# Set a float to define the conversion factor from Bohr radius to Angstrom.
# Recommendation: 0.529177210903 (CODATA 2018).
Bohr2Ang = 0.529177210903
# Set a float to define the conversion factor from Hartree to electronvolt.
# Recommendation: 27.211386245988 (CODATA 2018).
Hartree2eV = 27.211386245988

# The following does not required any modifications for regular usage.
# Set an array of floats which defines the tolerances in increasing order affecting the selection of
# extrapolated structures. The second entry specifies the threshold for the sum of the normalised
# symmetry function value extrapolations of the first selected extrapolated structure. If less than
# 0.1% of the simulations include such an extrapolation, the first entry is used instead. The
# following entries specify tested thresholds for the second selected structure. The initially used
# one is specified by initial_tolerance (initial_tolerance = 5 means sixth entry as Python starts
# counting at 0). Its value is increased if the normalised symmetry function value extrapolations of
# first and second selected structures overlap too much or if the minimum time step separation
# criterium is not fulfilled. The tolerance will be decreased if no large extrapolations are found
# or if the structure does not obey the geometrical rules specified above. In this way the entries
# between the second and by initial_tolerance specified entry of the array can be selected. The
# given values yielded good performance in all previous tests. If very small extrapolations are a
# problem, reduce the first two values. If there is a large gap between the small and large
# extrapolations, reduce the third to last values. You can also increase the number and density of
# the third to last entry to be more sensitive but with the drawback of a reduced performance. The
# value of initial_tolerance has to be higher than 1.
# Recommendation: [0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7,
# 0.8, 0.9, 1.0].
tolerances = [0.000001, 0.00001, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# Recommendation: 5.
initial_tolerance = 5

####################################################################################################

####################################################################################################

'''
Program start: No modifications are required for regular usage from here on.

'''

try:
  mode = sys.argv[1]
except IndexError:
  pass

exec(open(os.path.dirname(os.path.abspath(__file__))+'/lib/RuNNerActiveLearn_'+str(mode)+'.py',encoding="utf-8").read())

####################################################################################################

####################################################################################################
