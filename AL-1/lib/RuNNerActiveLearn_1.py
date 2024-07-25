
####################################################################################################

####################################################################################################

'''
RuNNerActiveLearn_1.py
written by Marco Eckhoff (2nd July 2021)

generates LAMMPS input files.

No modifications are required for regular usage.
'''

####################################################################################################

####################################################################################################

def check_input(integrators, pressures, N_steps, barostat_option, atom_style, dump_lammpstrj, min_timestep_separation_interpolation, element_types, masses, max_len_joblist, comment_name_keyword, structure_names, structure_selection, timestep):
  '''
  '''
  for integrator in integrators:
    if integrator!='nve' and integrator!='nvt' and integrator!='npt':
      print('ERROR: Integrator {0} is not implemented in RuNNerActiveLearn_1.py.'.format(integrator))
      exit()
  if not list(pressures) and 'npt' in integrators:
    print('ERROR: Integrator npt requires to specify at least one value for pressure.')
    exit()
  if barostat_option!='tri' and barostat_option!='aniso' and barostat_option!='iso':
    print('ERROR: Barostat option {0} is not implemented in RuNNerActiveLearn_1.py.'.format(barostat_option))
    exit()
  if atom_style!='atomic' and atom_style!='full':
    print('ERROR: Atom style {0} is not implemented in RuNNerActiveLearn_1.py.'.format(atom_style))
    exit()
  if N_steps%dump_lammpstrj!=0:
    print('ERROR: N_steps has to be a multiple of dump_lammpstrj ({0}!=N*{1}).'.format(N_steps, dump_lammpstrj))
    exit()
  if dump_lammpstrj<np.array(min_timestep_separation_interpolation).min():
    print('ERROR: The extrapolation free structures would be stored only every {0}th time step, but the minimum time step separation of interpolated structures is set to {1} time steps.'.format(dump_lammpstrj, np.array(min_timestep_separation_interpolation).min()))
    exit()
  if len(element_types)!=len(masses):
    print('ERROR: The number of given element types is not equal to the number of given masses ({0}!={1}).'.format(len(element_types), len(masses)))
    exit()
  if not max_len_joblist>=0:
    print('ERROR: The maximal length of the job list has to be set to 0 (which means infinity) or a positive integer number.'.format(max_len_joblist))
    exit()
  if (comment_name_keyword==None and structure_names!=None) or (comment_name_keyword!=None and structure_names==None):
    print('ERROR: If comment_name_keyword or structure_names is set to None the other one has to be set to None as well.')
    exit()
  if structure_names!=None:
    if list(structure_names):
      for structure_name in structure_names:
        if structure_name==None:
          print('ERROR: Individual structure names cannot be set to None. You have to specify an array of structure names or use structure_names = None.')
          exit()
    else:
      print('ERROR: structure_names has to be set to None or an array of structure names.')
      exit()
  if structure_names==None:
    if len(structure_selection)!=1:
      print('ERROR: As no structure names are given, exactly one setting for structure_selection is required.')
      exit()
    if structure_selection[0][0]<0 or structure_selection[0][1]<1:
      print('ERROR: The settings of structure_selection are not reasonable (every {0}th structure starting with the {0}th one).'.format(structure_selection[0][1], structure_selection[0][0]))
      exit()
  else:
    n_structure_names = len(structure_names)
    if len(structure_selection)==1:
      if structure_selection[0][0]<0 or structure_selection[0][1]<1:
        print('ERROR: The settings of structure_selection are not reasonable (every {0}th structure starting with the {0}th one).'.format(structure_selection[0][1], structure_selection[0][0]))
        exit()
      structure_selection = [deepcopy(structure_selection[0]) for i in range(n_structure_names)]
    elif len(structure_selection)!=n_structure_names:
      print('ERROR: Structure name dependent settings for structure_selection are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    else:
      for i in range(n_structure_names):
        if structure_selection[i][0]<0 or structure_selection[i][1]<1:
          print('ERROR: The settings of structure_selection are not reasonable (every {0}th structure starting with the {0}th one).'.format(structure_selection[i][1], structure_selection[i][0]))
          exit()
  if timestep>0.01:
    print('WARNING: Very large timestep of {0} ps.'.format(timestep))

  return structure_selection

####################################################################################################

def read_input_data(comment_name_keyword, comment_name_index, comment_name_separator, periodic, runner_cutoff):
  '''
  '''
  names = []
  lattices = []
  elements = []
  xyzs = []
  qs = []

  with open('RuNNer/input.data',encoding="utf-8") as f:
    for line in f.readlines():
      line = line.strip()
      if line.startswith('atom'):
        line = line.split()
        elements[-1].append(line[4])
        xyzs[-1].append([line[1], line[2], line[3]])
        qs[-1].append(line[5])
      elif line.startswith('lattice'):
        lattices[-1].append(line.split()[1:4])
      elif line.startswith('begin'):
        lattices.append([])
        elements.append([])
        xyzs.append([])
        qs.append([])
      elif line.startswith('end'):
        if not elements[-1]:
          print('ERROR: For some of the structures the definition of the atoms is incomplete or missing.')
          exit()
        xyzs[-1] = np.array(xyzs[-1]).astype(float)*Bohr2Ang
        qs[-1] = np.array(qs[-1]).astype(float)
        if periodic:
          if len(lattices[-1])==3:
            lattices[-1] = np.array(lattices[-1]).astype(float)*Bohr2Ang
          else:
            print('ERROR: The periodic keyword is set to True but for some of the structures the definition of the lattice is incomplete or missing.')
            exit()
        else:
          if lattices[-1]:
            print('ERROR: The periodic keyword is set to False but for some of the structures a definition of a lattice exists.')
            exit()
          else:
            lattices[-1] = np.array([[xyzs[-1][:,0].max()-xyzs[-1][:,0].min()+2*runner_cutoff*Bohr2Ang, 0.0, 0.0], [0.0, xyzs[-1][:,1].max()-xyzs[-1][:,1].min()+2*runner_cutoff*Bohr2Ang, 0.0], [0.0, 0.0, xyzs[-1][:,2].max()-xyzs[-1][:,2].min()+2*runner_cutoff*Bohr2Ang]])
      else:
        if comment_name_keyword!=None:
          if line.startswith(comment_name_keyword):
            names.append(line.split()[comment_name_index].split(comment_name_separator)[0])

  names = np.array(names)
  lattices = np.array(lattices)
  elements = np.array(elements)
  xyzs = np.array(xyzs)
  qs = np.array(qs)

  return names, lattices, elements, xyzs, qs

####################################################################################################

def write_input_lammps(path, seed, temperature, pressure, timestep, N_steps, integrator, barostat_option, runner_cutoff, element_types, atom_style, periodic):
  '''
  '''
  runner_cutoff = round(runner_cutoff*Bohr2Ang, 12)
  cflength = round(1.0/Bohr2Ang, 12)
  cfenergy = round(1.0/Hartree2eV, 15)
  elements_string = ''
  for element_type in element_types:
    elements_string += element_type+' '

  input_lammps = 'variable temperature equal {0}\n'.format(float(temperature))
  if integrator=='npt':
    input_lammps += 'variable pressure equal {0}\n'.format(float(pressure))
  input_lammps += 'variable N_steps equal {0}\n'.format(N_steps)\
                + 'variable seed equal {0}\n\n'.format(seed)
  input_lammps += 'units metal\n'\
                + 'boundary p p p\n'\
                + 'atom_style {0}\n'.format(atom_style)\
                + 'read_data structure.lammps\n'\
                + 'pair_style nnp dir RuNNer showew yes resetew no maxew 750 showewsum 0 cflength {0} cfenergy {1}\n'.format(cflength, cfenergy)\
                + 'pair_coeff * * {0}\n'.format(runner_cutoff)\
                + 'timestep {0}\n'.format(timestep)
  if integrator=='nve':
    input_lammps += 'fix int all nve\n'
  elif integrator=='nvt':
    input_lammps += 'fix int all nvt temp ${{temperature}} ${{temperature}} {0}\n'.format(timestep*100)
  elif integrator=='npt':
    input_lammps += 'fix int all npt temp ${{temperature}} ${{temperature}} {0} {1} ${{pressure}} ${{pressure}} {2} fixedpoint 0.0 0.0 0.0\n'.format(timestep*100, barostat_option, timestep*1000)
  input_lammps += 'thermo 1\n'\
                + 'variable thermo equal 0\n'\
                + 'thermo_style custom v_thermo step time temp epair etotal fmax fnorm press cella cellb cellc cellalpha cellbeta cellgamma density\n'\
                + 'thermo_modify format line "thermo %8d %10.4f %8.3f %15.5f %15.5f %9.4f %9.4f %9.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %8.5f"\n'
  if periodic:
    if atom_style=='atomic':
      input_lammps += 'dump lammpstrj all custom 1 structure.lammpstrj id element x y z\n'
    elif atom_style=='full':
      input_lammps += 'dump lammpstrj all custom 1 structure.lammpstrj id element x y z q\n'
    input_lammps += 'dump_modify lammpstrj pbc yes sort id element {0}\n'.format(elements_string[:-1])
  else:
    if atom_style=='atomic':
      input_lammps += 'dump lammpstrj all custom 1 structure.lammpstrj id element xu yu zu\n'
    elif atom_style=='full':
      input_lammps += 'dump lammpstrj all custom 1 structure.lammpstrj id element xu yu zu q\n'
    input_lammps += 'dump_modify lammpstrj pbc no sort id element {0}\n'.format(elements_string[:-1])
  input_lammps += 'velocity all create ${temperature} ${seed}\n\n'

  with open('mode1/'+path+'/input.lammps', 'w',encoding="utf-8") as f:
    f.write(input_lammps)

  subprocess.Popen('cat simulation.lammps >> mode1/'+path+'/input.lammps', shell=True)

####################################################################################################

def write_structure_lammps(path, lattice, element, xyz, q, element_types, masses, barostat_option, atom_style):
  '''
  '''
  lattice_lammps = transform_lattice(lattice)

  structure_lammps = 'RuNNerActiveLearn\n\n'\
                   + '{0} atoms\n\n'.format(len(element))\
                   + '{0} atom types\n\n'.format(len(element_types))\
                   + '0.0 {0} xlo xhi\n'.format(round(lattice_lammps[0], 5))\
                   + '0.0 {0} ylo yhi\n'.format(round(lattice_lammps[1], 5))\
                   + '0.0 {0} zlo zhi\n'.format(round(lattice_lammps[2], 5))
  if barostat_option=='tri':
    structure_lammps += '{0} {1} {2} xy xz yz\n'.format(round(lattice_lammps[3], 5), round(lattice_lammps[4],5), round(lattice_lammps[5], 5))
  structure_lammps += '\nMasses\n\n'
  for i in range(len(masses)):
    structure_lammps += '{0} {1}\n'.format(i+1, masses[i])
  structure_lammps += '\nAtoms\n\n'

  with open('mode1/'+path+'/structure.lammps', 'w',encoding="utf-8") as f:
    f.write(structure_lammps)
    if atom_style=='atomic':
      for i in range(len(element)):
        f.write('{0:4d} {1} {2:9.5f} {3:9.5f} {4:9.5f}\n'.format(i+1, element_types.index(element[i])+1, xyz[i][0], xyz[i][1], xyz[i][2]))
    elif atom_style=='full':
      for i in range(len(element)):
        f.write('{0:4d} 1 {1} {2:6.3f} {3:9.5f} {4:9.5f} {5:9.5f}\n'.format(i+1, element_types.index(element[i])+1, round(q[i], 3), round(xyz[i][0], 5), round(xyz[i][1], 5), round(xyz[i][2], 5)))

####################################################################################################

def transform_lattice(lattice):
  '''
  '''
  a = np.linalg.norm(lattice[0])
  b = np.linalg.norm(lattice[1])
  c = np.linalg.norm(lattice[2])
  cos_alpha = np.dot(lattice[1], lattice[2])/b/c
  cos_beta = np.dot(lattice[0], lattice[2])/a/c
  cos_gamma = np.dot(lattice[0], lattice[1])/a/b
  xy = b*cos_gamma
  lx = a
  xz = c*cos_beta
  ly = np.sqrt(b**2-xy**2)
  yz = (b*c*cos_alpha-xy*xz)/ly
  lz = np.sqrt(c**2-xz**2-yz**2)

  return [lx, ly, lz, xy, xz, yz]

####################################################################################################

####################################################################################################

structure_selection = check_input(integrators, pressures, N_steps, barostat_option, atom_style, dump_lammpstrj, min_timestep_separation_interpolation, element_types, masses, max_len_joblist, comment_name_keyword, structure_names, structure_selection, timestep)
if os.path.isdir('mode1'):
  print('ERROR: Path mode1 already exists. Please remove old directory first if you would like to recreate it.')
  exit()
subprocess.Popen('mkdir mode1', shell=True).wait()

names_all, lattices_all, elements_all, xyzs_all, qs_all = read_input_data(comment_name_keyword, comment_name_index, comment_name_separator, periodic, runner_cutoff)

if max_len_joblist==0:
  joblist_name = 'joblist_mode1.dat'
  with open(joblist_name, 'w',encoding="utf-8") as f:
    f.write('')
if structure_names==None:
  structure_names = [None]
pressures_npt = pressures
n_simulations = 0
n_previous_simulations = 0
counter = 0

for i in range(len(structure_names)):
  if structure_names[i]==None:
    names = names_all
    lattices = lattices_all
    elements = elements_all
    xyzs = xyzs_all
    qs = qs_all
  else:
    names = names_all[names_all==structure_names[i]]
    lattices = lattices_all[names_all==structure_names[i]]
    elements = elements_all[names_all==structure_names[i]]
    xyzs = xyzs_all[names_all==structure_names[i]]
    qs = qs_all[names_all==structure_names[i]]
    print('Structure name: {0}'.format(structure_names[i]))
  structure_selection[i][0] = structure_selection[i][0]%structure_selection[i][1]
  print('Starting from the {0}th structure every {1}th structure of the input.data file is used.'.format(structure_selection[i][0], structure_selection[i][1]))
  n_structures = len(lattices)
  n_npt = int(np.array([1 for j in integrators if j=='npt']).sum())
  repeatitions = max(1, int(float(n_structures)/2/len(integrators)/len(temperatures)/((len(pressures)-1)*n_npt/len(integrators)+1)/structure_selection[i][1]))
  print('The given variations of the settings are repeated {0} times'.format(repeatitions))

  for x in range(repeatitions):
    for HDNNP in ['1', '2']:
      for integrator in integrators:
        for temperature in temperatures:
          if integrator!='npt':
            pressures = [0]
          else:
            pressures = pressures_npt
          for pressure in pressures:
            if n_structures//structure_selection[i][1]<=counter:
              n_simulations += counter
              counter = 0
              print('WARNING: The structures of the input.data file are used more than once.')
              if structure_selection[i][1]>1:
                structure_selection[i][0] = (structure_selection[i][0]+1)%structure_selection[i][1]
                print('Try to avoid this by start from the {0}th structure and using again every {1}th structure.'.format(structure_selection[i][0], structure_selection[i][1]))
            selection = counter*structure_selection[i][1]+structure_selection[i][0]
            if comment_name_keyword!=None:
              if integrator=='npt':
                path = names[selection]+'_'+integrator+'_hdnnp'+HDNNP+'_t'+str(temperature)+'_p'+str(pressure)+'_'+str(seed)
              else:
                path = names[selection]+'_'+integrator+'_hdnnp'+HDNNP+'_t'+str(temperature)+'_'+str(seed)
            else:
              path = integrator+'_hdnnp'+HDNNP+'_t'+str(temperature)+'_p'+str(pressure)+'_'+str(seed)
            if os.path.isdir('mode1/'+path):
              print('ERROR: Path mode1/{0} already exists. Please remove old directories first if you would like to recreate them.'.format(path))
              exit()
            subprocess.Popen('mkdir mode1/'+path, shell=True).wait()
            write_input_lammps(path, seed, temperature, pressure, timestep, N_steps, integrator, barostat_option, runner_cutoff, element_types, atom_style, periodic)
            write_structure_lammps(path, lattices[selection], elements[selection], xyzs[selection], qs[selection], element_types, masses, barostat_option, atom_style)
            subprocess.Popen('mkdir mode1/'+path+'/RuNNer', shell=True).wait()
            subprocess.Popen('cp -i RuNNer/HDNNP_'+HDNNP+'/* mode1/'+path+'/RuNNer', shell=True)
            if max_len_joblist!=0 and (n_simulations+counter)%max_len_joblist==0:
              joblist_name = 'joblist_mode1_'+str((n_simulations+counter)//max_len_joblist+1)+'.dat'
              with open(joblist_name, 'w',encoding="utf-8") as f:
                f.write('')
            with open(joblist_name, 'a',encoding="utf-8") as f:
              f.write('{0}\n'.format(path))
            seed += 1
            counter += 1

  if structure_names[i]!=None:
    n_simulations += counter
    counter = 0
    print('Input was generated for {0} simulations.'.format(n_simulations-n_previous_simulations))
    n_previous_simulations = n_simulations

if structure_names[0]==None:
  n_simulations += counter
  print('Input was generated for {0} simulations.'.format(n_simulations))

####################################################################################################

####################################################################################################
