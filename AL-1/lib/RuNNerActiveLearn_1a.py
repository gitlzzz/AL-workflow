####################################################################################################

####################################################################################################

'''
RuNNerActiveLearn_1a.py
written by Marco Eckhoff (2nd July 2021)

preanalyses the extrapolations and reduces the size of the structure.lammpstrj file.

No modifications are required for regular usage.
'''

####################################################################################################

####################################################################################################

def read_log(dump_lammpstrj):
  '''
  '''
  with open('log.lammps',encoding="utf-8") as f:
    data = [line for line in f.readlines()]
  counter = 0
  n_lines = len(data)
  while counter<n_lines and not data[counter].startswith('**********'):
    counter += 1

  extrapolation = False
  i = counter
  while i<n_lines and not data[i].startswith('### NNP EXTRAPOLATION WARNING ###'):
    i += 1
  if i<n_lines:
    extrapolation = True
  i -= 1
  while i>counter and not data[i].startswith('thermo'):
    i -= 1
  if extrapolation:
    extrapolation_free_lines = i
    if i>counter:
      extrapolation_free_timesteps = int(data[i].split()[1])
    else:
      extrapolation_free_timesteps = -1
  else:
    extrapolation_free_lines = -1
    extrapolation_free_timesteps = int(data[i].split()[1])

  data = [int(line.split()[1]) if line.startswith('thermo') else -1 for line in data[counter:] if line.startswith('thermo') or line.startswith('### NNP EXTRAPOLATION WARNING ###')]
  timesteps = np.unique(np.array([data[i] for i in range(1,len(data)) if data[i]!=-1 and (data[i]%dump_lammpstrj==0 or data[i-1]==-1)]))

  return timesteps, extrapolation_free_lines, extrapolation_free_timesteps

####################################################################################################

def read_lammpstrj(timesteps):
  '''
  '''
  structures = []
  i = 0
  n_timesteps = len(timesteps)
  with open('structure.lammpstrj',encoding="utf-8") as f:
    line = f.readline()
    while line and i<n_timesteps:
      while not line.startswith('ITEM: TIMESTEP') and line:
        line = f.readline()
      line = f.readline()
      if timesteps[i]==int(line.strip()):
        structures.append('ITEM: TIMESTEP\n')
        while not line.startswith('ITEM: TIMESTEP') and line:
          structures.append(line)
          line = f.readline()
        i += 1

  i = 1
  n_lines = len(structures)
  while i<n_lines and not structures[i].startswith('ITEM: TIMESTEP'):
    i += 1
  structure_lines = i 

  return structures, structure_lines

####################################################################################################

def write_lammpstrj(structures):
  '''
  '''
  with open('structure.lammpstrj', 'w',encoding="utf-8") as f:
    for line in structures:
      f.write(line)

####################################################################################################

def write_extrapolation(extrapolation_free_timesteps, extrapolation_free_lines, dump_lammpstrj, structure_lines, last_timestep):
  '''
  '''
  with open('extrapolation.dat', 'w',encoding="utf-8") as f:
    f.write('extrapolation_free_initial_time_steps: {0}\nlines_before_first_extrapolation: {1}\ntimesteps_between_non_extrapolated_structures: {2}\nlines_per_structure: {3}\nlast_timestep: {4}'.format(extrapolation_free_timesteps, extrapolation_free_lines, dump_lammpstrj, structure_lines, last_timestep))

####################################################################################################

####################################################################################################

timesteps, extrapolation_free_lines, extrapolation_free_timesteps = read_log(dump_lammpstrj)

structures, structure_lines = read_lammpstrj(timesteps)

write_lammpstrj(structures)

write_extrapolation(extrapolation_free_timesteps, extrapolation_free_lines, dump_lammpstrj, structure_lines, timesteps[-1])

####################################################################################################

####################################################################################################
