###################################################################################################

'''
RuNNerActiveLearn_2.py
written by Marco Eckhoff (2nd July 2021)
modify to 2.1.4 by Zan 30 Apr 2022
selects candidates for new reference structures from simulation outputs, recalculates the selected
structures employing the two HDNNPs, and compares their predictions of energies and forces to find
missing structures which exhibit energy and force errors above the given thresholds. Statistics of
the extrapolated structures are written as well.

No modifications are required for regular usage.
'''

####################################################################################################

####################################################################################################

def check_input(element_types, masses, structure_names, d_mins, min_timestep_separation_extrapolation, timestep_separation_interpolation_checks, min_timestep_separation_interpolation, delta_E, delta_F, all_extrapolated_structures, max_extrapolated_structures, max_interpolated_structures_per_simulation, exceptions, tolerances, initial_tolerance):
  '''
  '''
  if len(element_types)!=len(masses):
    print('ERROR: The number of given element types is not equal to the number of given masses ({0}!={1}).'.format(len(element_types), len(masses)))
    exit()
  if structure_names==None:
    if not (1==len(d_mins)==len(min_timestep_separation_extrapolation)==len(timestep_separation_interpolation_checks)==len(min_timestep_separation_interpolation)==len(delta_E)==len(delta_F)==len(all_extrapolated_structures)==len(exceptions)):
      print('ERROR: As no structure names are given, one setting for each structure name dependent setting is required.')
      exit()
  else:
    if not list(structure_names):
      print('ERROR: structure_names has to be set to None or an array of structure names.')
      exit()
    n_structure_names = len(structure_names)
    if len(d_mins)==1:
      d_mins *= n_structure_names
    elif len(d_mins)!=n_structure_names:
      print('ERROR: Structure name dependent settings for d_mins are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(min_timestep_separation_extrapolation)==1:
      min_timestep_separation_extrapolation *= n_structure_names
    elif len(min_timestep_separation_extrapolation)!=n_structure_names:
      print('ERROR: Structure name dependent settings for min_timestep_separation_extrapolation are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(timestep_separation_interpolation_checks)==1:
      timestep_separation_interpolation_checks *= n_structure_names
    elif len(timestep_separation_interpolation_checks)!=n_structure_names:
      print('ERROR: Structure name dependent settings for timestep_separation_interpolation_checks are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(min_timestep_separation_interpolation)==1:
      min_timestep_separation_interpolation *= n_structure_names
    elif len(min_timestep_separation_interpolation)!=n_structure_names:
      print('ERROR: Structure name dependent settings for min_timestep_separation_interpolation are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(delta_E)==1:
      delta_E *= n_structure_names
    elif len(delta_E)!=n_structure_names:
      print('ERROR: Structure name dependent settings for delta_E are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(delta_F)==1:
      delta_F *= n_structure_names
    elif len(delta_F)!=n_structure_names:
      print('ERROR: Structure name dependent settings for delta_F are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(all_extrapolated_structures)==1:
      all_extrapolated_structures *= n_structure_names
    elif len(all_extrapolated_structures)!=n_structure_names:
      print('ERROR: Structure name dependent settings for all_extrapolated_structures are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(max_extrapolated_structures)==1:
      max_extrapolated_structures *= n_structure_names
    elif len(max_extrapolated_structures)!=n_structure_names:
      print('ERROR: Structure name dependent settings for max_extrapolated_structures are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(max_interpolated_structures_per_simulation)==1:
      max_interpolated_structures_per_simulation *= n_structure_names
    elif len(max_interpolated_structures_per_simulation)!=n_structure_names:
      print('ERROR: Structure name dependent settings for max_interpolated_structures_per_simulation are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if len(exceptions)==1:
      exceptions *= n_structure_names
    elif len(exceptions)!=n_structure_names:
      print('ERROR: Structure name dependent settings for exceptions are not given for every structure name or there are too many settings for the given structure names. Also there is not given one value which could be used for all structures names.')
      exit()
    if not (np.array(max_extrapolated_structures)>=0).all():
      print('ERROR: The value of max_extrapolated_structures has to be an integer equal or higher than 0.')
      exit()
  if (np.array(timestep_separation_interpolation_checks)*5>=N_steps).any():
    print('ERROR: The time step separation between two interpolation checks has to be smaller than a fifth of the number of MD steps.')
    exit()
  if (np.array(timestep_separation_interpolation_checks)<np.array(min_timestep_separation_interpolation)).any():
    print('ERROR: The time step separation between two interpolation checks is set to a smaller value than the minimal time step separation between two interpolated structures.')
    exit()
  if initial_tolerance<=1:
    print('ERROR: The value of initial_tolerance has to be higher than 1.')
    exit()
  if len(tolerances)<=initial_tolerance:
    print('ERROR: There are not enough tolerance values as initial_tolerance results in an index error.')
    exit()

  return d_mins, min_timestep_separation_extrapolation, timestep_separation_interpolation_checks, min_timestep_separation_interpolation, delta_E, delta_F, all_extrapolated_structures, max_extrapolated_structures, max_interpolated_structures_per_simulation, exceptions

####################################################################################################

def get_paths(structure_name):
  '''
  '''
  if structure_name!='':
    try:
      paths = str(subprocess.check_output('ls mode1 | grep '+structure_name+'_ | grep -e _nve_hdnnp -e _nvt_hdnnp -e _npt_hdnnp', stderr=subprocess.STDOUT, shell=True).decode()).strip().split('\n')
    except subprocess.CalledProcessError:
      print('ERROR: Simulations with the structure name {0} were not found.'.format(structure_name))
      exit()
  else:
    paths = str(subprocess.check_output('ls mode1 | grep -e nve_hdnnp -e nvt_hdnnp -e npt_hdnnp', stderr=subprocess.STDOUT, shell=True).decode()).strip().split('\n')

  finished = []
  for i in range(len(paths)):
    if os.path.isfile('mode1/'+paths[i]+'/extrapolation.dat'):
      finished.append(i)
    else:
      print('Simulation {0} is not finished.'.format(paths[i]))
  paths = np.array(paths)[finished]
  if len(paths)==0:
    print('ERROR: None of the {0} simulations finished.'.format(structure_name))
    exit()

  return paths

####################################################################################################

def read_extrapolation(path):
  '''
  '''
  with open('mode1/'+path+'/extrapolation.dat',encoding="utf-8") as f:
    extrapolation_data = np.array([line.strip().split() for line in f.readlines()])[:,1].astype(int)

  return extrapolation_data

####################################################################################################

def read_log_format(path):
  '''
  '''
  with open('mode1/'+path+'/log.lammps',encoding="utf-8") as f:
    data = [line for line in f.readlines()]
  counter = 0
  n_lines = len(data)
  while counter<n_lines and not data[counter].startswith('**********'):
    counter += 1
  if counter<n_lines:
    if data[counter+2].startswith('   NNP LIBRARY v2.0.0'):
      extrapolation_format = 'v2.0.0'
    elif data[counter+5].startswith('n²p² version  (from git): v2.1.4'):
      extrapolation_format = 'v2.1.4'
    else:
      print('ERROR: n2p2 extrapolation warning format cannot be identified in the file {0}/log.lammps. Known formats are corresponding to n2p2 v2.0.0 and v2.1.4.'.format(path))
      exit()
  else:
    print('ERROR: n2p2 extrapolation warning format cannot be identified in the file {0}/log.lammps. Known formats are corresponding to n2p2 v2.0.0 and v2.1.4.'.format(path))
    exit()

  return extrapolation_format

####################################################################################################

def read_log(path, extrapolation_data, tolerances, extrapolation_format):
  '''
  '''
  if extrapolation_data[1]!=-1:
    with open('mode1/'+path+'/log.lammps',encoding="utf-8") as f:
      data = [line.strip() for line in f.readlines()][extrapolation_data[1]:-1]
    if extrapolation_format=='v2.0.0':
      data = np.array([[float(line.split()[1])]+[np.nan,np.nan,np.nan,np.nan,np.nan] if line.startswith('thermo') else [np.nan]+list(np.array(line.split())[[12,14,16,8,10]].astype(float)) for line in data if line.startswith('thermo') or line.startswith('### NNP')])
    elif extrapolation_format=='v2.1.4':
      data = np.array([[float(line.split()[1])]+[np.nan,np.nan,np.nan,np.nan,np.nan] if line.startswith('thermo') else [np.nan]+list(np.array(line.split())[[16,18,20,8,12]].astype(float)) for line in data if line.startswith('thermo') or line.startswith('### NNP')])
    if np.isnan(data[-1,0]):
      data = data[:-np.argmax(np.isfinite(data[:,0][::-1]))]
    if np.isnan(data[0,0]):
      print('WARNING: Extrapolation occurred already in the first time step in {0}.'.format(path))
      data = np.concatenate((np.array([[-1.0]+[np.nan,np.nan,np.nan,np.nan,np.nan]]),data), axis=0)
    extrapolation = np.absolute((data[:,1]-data[:,2])/(data[:,3]-data[:,2])-0.5)-0.5
    for i in range(1,len(extrapolation)):
      if np.isfinite(extrapolation[i]) and np.isfinite(extrapolation[i-1]):
        extrapolation[i] += extrapolation[i-1]
    extrapolation_indices = []
    with warnings.catch_warnings():
      warnings.simplefilter('ignore', category=RuntimeWarning)
      for tolerance in tolerances:
        extrapolation_indices.append(np.argmax(extrapolation>tolerance))
    extrapolation_timestep = []
    extrapolation_value = []
    extrapolation_statistic = []
    for i in range(len(tolerances)):
      if extrapolation_indices[i]>0:
        j = 1
        while np.isnan(data[extrapolation_indices[i]+j,0]):
          j += 1
        extrapolation_timestep.append(int(data[extrapolation_indices[i]+j,0]))
        extrapolation_value.append(extrapolation[extrapolation_indices[i]+j-1])
        extrapolation_statistic.append([])
        j -= 1
        extrapolation_statistic[-1].append(data[extrapolation_indices[i]+j,[1,4,5]])
        extrapolation_statistic[-1][-1][0] = (data[extrapolation_indices[i]+j,1]-data[extrapolation_indices[i]+j,2])/(data[extrapolation_indices[i]+j,3]-data[extrapolation_indices[i]+j,2])-0.5
        if extrapolation_statistic[-1][-1][0]<0:
          extrapolation_statistic[-1][-1][0] += 0.5
        else:
          extrapolation_statistic[-1][-1][0] -= 0.5
        j -= 1
        while np.isnan(data[extrapolation_indices[i]+j,0]):
          extrapolation_statistic[-1].append(data[extrapolation_indices[i]+j,[1,4,5]])
          extrapolation_statistic[-1][-1][0] = (data[extrapolation_indices[i]+j,1]-data[extrapolation_indices[i]+j,2])/(data[extrapolation_indices[i]+j,3]-data[extrapolation_indices[i]+j,2])-0.5
          if extrapolation_statistic[-1][-1][0]<0:
            extrapolation_statistic[-1][-1][0] += 0.5
          else:
            extrapolation_statistic[-1][-1][0] -= 0.5
          j -= 1
      else:
        extrapolation_timestep.append(-1)
        extrapolation_value.append(0)
        extrapolation_statistic.append([None])
  else:
    extrapolation_timestep = len(tolerances)*[-1]
    extrapolation_value = len(tolerances)*[0]
    extrapolation_statistic = len(tolerances)*[None]

  return extrapolation_timestep, extrapolation_value, extrapolation_statistic

####################################################################################################

def get_timesteps(extrapolation_timesteps, extrapolation_values, extrapolation_data, min_timestep_separation_extrapolation, timestep_separation_interpolation_checks, min_timestep_separation_interpolation, tolerances, initial_tolerance):
  '''
  '''
  min_fraction = 0.001
  n_tolerances = len(tolerances)
  n_small = 2
  small = n_small-1
  if len(extrapolation_timesteps[:,small][extrapolation_timesteps[:,small]>=min_timestep_separation_extrapolation])<min_fraction*len(extrapolation_timesteps[:,small]):
    print('Only less than 0.1% of the simulations show an extrapolation if a tolerance of {0} is employed (the initial {1} time steps are neglected). The tolerance value is reduced to {2}.'.format(tolerances[small], min_timestep_separation_extrapolation, tolerances[small-1]))
    small -= 1
  if not (extrapolation_timesteps[:,small][extrapolation_timesteps[:,small]>=0]).any():
    print('There are no small extrapolations.')
    tolerance_indices = np.array(len(extrapolation_timesteps)*[-1])
  else:
    n_simulations_extrapolation = len(extrapolation_timesteps[:,small][extrapolation_timesteps[:,small]>=0])
    n_simulations = len(extrapolation_timesteps[:,small])
    print('Small extrapolations are present in {0} of {1} simulations ({2}%).'.format(n_simulations_extrapolation, n_simulations, round(100.0*n_simulations_extrapolation/n_simulations,2)))
    extrapolation_values_reduced = extrapolation_values[extrapolation_timesteps[:,small]==-1]
    mean_small = np.mean(extrapolation_values_reduced[:,small])
    std_small = np.std(extrapolation_values_reduced[:,small])
    criterium = mean_small+std_small+max(0,tolerances[small]-mean_small+std_small)
    while criterium>tolerances[initial_tolerance] and initial_tolerance<n_tolerances:
      initial_tolerance += 1
    while not (extrapolation_timesteps[:,initial_tolerance][extrapolation_timesteps[:,initial_tolerance]>=min_timestep_separation_extrapolation]).any() and initial_tolerance>n_small:
      print('There are no large extrapolations for a tolerance of {0} (the initial {1} time steps are neglected). The tolerance value is reduced to {2}.'.format(tolerances[initial_tolerance], min_timestep_separation_extrapolation, tolerances[initial_tolerance-1]))
      initial_tolerance -= 1
    if initial_tolerance==n_small:
      if not (extrapolation_timesteps[:,initial_tolerance][extrapolation_timesteps[:,initial_tolerance]>=0]).any():
        print('There are no large extrapolations.')
    extra_steps = (extrapolation_timesteps[:,n_small:].T-extrapolation_timesteps[:,small]).T
    extra_steps[extra_steps<0] = min_timestep_separation_extrapolation+1
    extra_steps_reduced = extra_steps[extrapolation_timesteps[:,small]!=-1]
    tolerance_indices = initial_tolerance*np.ones(len(extra_steps), dtype=int)
    tolerance_indices[extrapolation_timesteps[:,small]==-1] = -1
    tolerance_indices_reduced = tolerance_indices[extrapolation_timesteps[:,small]!=-1]
    for i in range(initial_tolerance-n_small, n_tolerances-n_small):
      tolerance_indices_reduced[extra_steps_reduced[:,i]<min_timestep_separation_extrapolation] += 1
    tolerance_indices_reduced[tolerance_indices_reduced>=n_tolerances] = -1
    tolerance_indices[extrapolation_timesteps[:,small]!=-1] = tolerance_indices_reduced
    tolerance_indices[tolerance_indices>=small] -= n_small

  selected_timesteps = []
  smalls = small*np.ones(len(extrapolation_timesteps), dtype=int)
  min_interpolated_structure_checks = 3
  for i in range(len(extrapolation_timesteps)):
    if extrapolation_timesteps[i][small]<0 and extrapolation_data[i][4]!=N_steps:
      print('WARNING: A simulation ended due to too many extrapolations but no one of these was larger than the tolerance of {0}. If this message is printed several times you should consider to reduce the first and second entry of tolerances.'.format(tolerances[small]))
      if small>0 and extrapolation_timesteps[i][small-1]>=min_timestep_separation_extrapolation:
        smalls[i] = small-1
        print('With the reduction of the tolerance to {0} an extrapolated structure could be found in this case.'.format(tolerances[smalls[i]]))
    if extrapolation_timesteps[i][smalls[i]]>=0:
      if extrapolation_timesteps[i][smalls[i]]>(min_interpolated_structure_checks+2)*timestep_separation_interpolation_checks:
        selected_timesteps.append(list(range(2*timestep_separation_interpolation_checks, extrapolation_timesteps[i][smalls[i]]-timestep_separation_interpolation_checks+1, timestep_separation_interpolation_checks))+[extrapolation_timesteps[i][smalls[i]], deepcopy(extrapolation_timesteps[i][n_small:])])
      else:
        small_timestep_separation_interpolation_checks = ((extrapolation_timesteps[i][smalls[i]]//(min_interpolated_structure_checks+2))//extrapolation_data[i][2])*extrapolation_data[i][2]
        n_interpolation_checks = min_interpolated_structure_checks
        while small_timestep_separation_interpolation_checks<min_timestep_separation_interpolation and n_interpolation_checks>1:
          n_interpolation_checks -= 1
          small_timestep_separation_interpolation_checks = (extrapolation_timesteps[i][smalls[i]]//(n_interpolation_checks+2)//extrapolation_data[i][2])*extrapolation_data[i][2]
        if small_timestep_separation_interpolation_checks>min_timestep_separation_interpolation:
          selected_timesteps.append([j*small_timestep_separation_interpolation_checks for j in range(2,n_interpolation_checks+2)]+[extrapolation_timesteps[i][smalls[i]], deepcopy(extrapolation_timesteps[i][n_small:])])
        else:
          selected_timesteps.append([extrapolation_timesteps[i][smalls[i]], deepcopy(extrapolation_timesteps[i][n_small:])])
    else:
      if extrapolation_data[i][4]>(min_interpolated_structure_checks+2)*timestep_separation_interpolation_checks:
        selected_timesteps.append(list(range(2*timestep_separation_interpolation_checks, extrapolation_data[i][4]+1, timestep_separation_interpolation_checks))+[-1, (n_tolerances-n_small)*[-1]])
      else:
        small_timestep_separation_interpolation_checks = ((extrapolation_data[i][4]//(min_interpolated_structure_checks+2))//extrapolation_data[i][2])*extrapolation_data[i][2]
        n_interpolation_checks = min_interpolated_structure_checks
        while small_timestep_separation_interpolation_checks<min_timestep_separation_interpolation and n_interpolation_checks>1:
          n_interpolation_checks -= 1
          small_timestep_separation_interpolation_checks = (extrapolation_data[i][4]//(n_interpolation_checks+2)//extrapolation_data[i][2])*extrapolation_data[i][2]
        if small_timestep_separation_interpolation_checks>min_timestep_separation_interpolation:
          selected_timesteps.append([j*small_timestep_separation_interpolation_checks for j in range(2,n_interpolation_checks+2)]+[(extrapolation_data[i][4]//extrapolation_data[i][2])*extrapolation_data[i][2], -1, (n_tolerances-n_small)*[-1]])
          print('Included the last regularly dumped structure of the simulation as it ended due to too many extrapolations.')
        else:
          if (extrapolation_data[i][4]//extrapolation_data[i][2])*extrapolation_data[i][2]>=min_timestep_separation_extrapolation:
            selected_timesteps.append([(extrapolation_data[i][4]//extrapolation_data[i][2])*extrapolation_data[i][2], -1, (n_tolerances-n_small)*[-1]])
            print('Included the last regularly dumped structure of the simulation as it ended due to too many extrapolations.')
          else:
            selected_timesteps.append([-1, (n_tolerances-n_small)*[-1]])

  return selected_timesteps, tolerance_indices, smalls, n_small

#####################################################################################################

def get_structure(data):
  '''
  '''
  lat = np.array([data[5].split(), data[6].split(), data[7].split()]).astype(float)
  if data[4].startswith('ITEM: BOX BOUNDS xy xz yz pp pp pp'):
    lx = lat[0][1]-lat[0][0]-np.array([0.0,lat[0][2],lat[1][2],lat[0][2]+lat[1][2]]).max()+np.array([0.0,lat[0][2],lat[1][2],lat[0][2]+lat[1][2]]).min()
    ly = lat[1][1]-lat[1][0]-np.array([0.0,lat[2][2]]).max()+np.array([0.0,lat[2][2]]).min()
    lz = lat[2][1]-lat[2][0]
    lattice = [[lx, 0.0, 0.0], [lat[0][2], ly, 0.0], [lat[1][2], lat[2][2], lz]]
  else:
    lattice = [[(lat[0][1]-lat[0][0]), 0.0, 0.0], [0.0, (lat[1][1]-lat[1][0]), 0.0], [0.0, 0.0, (lat[2][1]-lat[2][0])]]

  atom_style = 'atomic'
  if data[8].startswith('ITEM: ATOMS id element x y z q') or data[8].startswith('ITEM: ATOMS id element xu yu zu q'):
    atom_style = 'full'
  data = np.array([line.split() for line in data[9:]])
  element = deepcopy(data[:,1])
  position = deepcopy(data[:,2:5]).astype(float)
  if atom_style=='full':
    charge = deepcopy(data[:,5]).astype(float)
  else:
    charge = np.zeros(len(element))

  return lattice, element, position, charge

####################################################################################################

def check_nearest_neighbours(lat, pos_i, pos_j, ii, d_min):
  '''
  '''
  if len(pos_i)==0 or len(pos_j)==0:
    return True, -1

  if pos_i.ndim==1:
    pos_i = np.array([pos_i])
  if pos_j.ndim==1:
    pos_j = np.array([pos_j])

  pos = np.array(deepcopy(pos_j))
  pos = np.concatenate((pos,
                        np.dstack((pos[:,0]-lat[0][0]-lat[1][0]-lat[2][0], pos[:,1]-lat[0][1]-lat[1][1]-lat[2][1], pos[:,2]-lat[0][2]-lat[1][2]-lat[2][2]))[0],
                        np.dstack((pos[:,0]-lat[0][0]-lat[1][0]          , pos[:,1]-lat[0][1]-lat[1][1]          , pos[:,2]-lat[0][2]-lat[1][2]          ))[0],
                        np.dstack((pos[:,0]-lat[0][0]-lat[1][0]+lat[2][0], pos[:,1]-lat[0][1]-lat[1][1]+lat[2][1], pos[:,2]-lat[0][2]-lat[1][2]+lat[2][2]))[0],
                        np.dstack((pos[:,0]-lat[0][0]          -lat[2][0], pos[:,1]-lat[0][1]          -lat[2][1], pos[:,2]-lat[0][2]          -lat[2][2]))[0],
                        np.dstack((pos[:,0]-lat[0][0]                    , pos[:,1]-lat[0][1]                    , pos[:,2]-lat[0][2]                    ))[0],
                        np.dstack((pos[:,0]-lat[0][0]          +lat[2][0], pos[:,1]-lat[0][1]          +lat[2][1], pos[:,2]-lat[0][2]          +lat[2][2]))[0],
                        np.dstack((pos[:,0]-lat[0][0]+lat[1][0]-lat[2][0], pos[:,1]-lat[0][1]+lat[1][1]-lat[2][1], pos[:,2]-lat[0][2]+lat[1][2]-lat[2][2]))[0],
                        np.dstack((pos[:,0]-lat[0][0]+lat[1][0]          , pos[:,1]-lat[0][1]+lat[1][1]          , pos[:,2]-lat[0][2]+lat[1][2]          ))[0],
                        np.dstack((pos[:,0]-lat[0][0]+lat[1][0]+lat[2][0], pos[:,1]-lat[0][1]+lat[1][1]+lat[2][1], pos[:,2]-lat[0][2]+lat[1][2]+lat[2][2]))[0],
                        np.dstack((pos[:,0]          -lat[1][0]-lat[2][0], pos[:,1]          -lat[1][1]-lat[2][1], pos[:,2]          -lat[1][2]-lat[2][2]))[0],
                        np.dstack((pos[:,0]          -lat[1][0]          , pos[:,1]          -lat[1][1]          , pos[:,2]          -lat[1][2]          ))[0],
                        np.dstack((pos[:,0]          -lat[1][0]+lat[2][0], pos[:,1]          -lat[1][1]+lat[2][1], pos[:,2]          -lat[1][2]+lat[2][2]))[0],
                        np.dstack((pos[:,0]                    -lat[2][0], pos[:,1]                    -lat[2][1], pos[:,2]                    -lat[2][2]))[0],
                        np.dstack((pos[:,0]                    +lat[2][0], pos[:,1]                    +lat[2][1], pos[:,2]                    +lat[2][2]))[0],
                        np.dstack((pos[:,0]          +lat[1][0]-lat[2][0], pos[:,1]          +lat[1][1]-lat[2][1], pos[:,2]          +lat[1][2]-lat[2][2]))[0],
                        np.dstack((pos[:,0]          +lat[1][0]          , pos[:,1]          +lat[1][1]          , pos[:,2]          +lat[1][2]          ))[0],
                        np.dstack((pos[:,0]          +lat[1][0]+lat[2][0], pos[:,1]          +lat[1][1]+lat[2][1], pos[:,2]          +lat[1][2]+lat[2][2]))[0],
                        np.dstack((pos[:,0]+lat[0][0]-lat[1][0]-lat[2][0], pos[:,1]+lat[0][1]-lat[1][1]-lat[2][1], pos[:,2]+lat[0][2]-lat[1][2]-lat[2][2]))[0],
                        np.dstack((pos[:,0]+lat[0][0]-lat[1][0]          , pos[:,1]+lat[0][1]-lat[1][1]          , pos[:,2]+lat[0][2]-lat[1][2]          ))[0],
                        np.dstack((pos[:,0]+lat[0][0]-lat[1][0]+lat[2][0], pos[:,1]+lat[0][1]-lat[1][1]+lat[2][1], pos[:,2]+lat[0][2]-lat[1][2]+lat[2][2]))[0],
                        np.dstack((pos[:,0]+lat[0][0]          -lat[2][0], pos[:,1]+lat[0][1]          -lat[2][1], pos[:,2]+lat[0][2]          -lat[2][2]))[0],
                        np.dstack((pos[:,0]+lat[0][0]                    , pos[:,1]+lat[0][1]                    , pos[:,2]+lat[0][2]                    ))[0],
                        np.dstack((pos[:,0]+lat[0][0]          +lat[2][0], pos[:,1]+lat[0][1]          +lat[2][1], pos[:,2]+lat[0][2]          +lat[2][2]))[0],
                        np.dstack((pos[:,0]+lat[0][0]+lat[1][0]-lat[2][0], pos[:,1]+lat[0][1]+lat[1][1]-lat[2][1], pos[:,2]+lat[0][2]+lat[1][2]-lat[2][2]))[0],
                        np.dstack((pos[:,0]+lat[0][0]+lat[1][0]          , pos[:,1]+lat[0][1]+lat[1][1]          , pos[:,2]+lat[0][2]+lat[1][2]          ))[0],
                        np.dstack((pos[:,0]+lat[0][0]+lat[1][0]+lat[2][0], pos[:,1]+lat[0][1]+lat[1][1]+lat[2][1], pos[:,2]+lat[0][2]+lat[1][2]+lat[2][2]))[0]), axis=0)

  if ii:
    select = 1
  else:
    select = 0

  for p in pos_i:
    d = np.dstack((pos[:,0]-p[0], pos[:,1]-p[1], pos[:,2]-p[2]))[0]
    d = np.sqrt(d[:,0]**2+d[:,1]**2+d[:,2]**2)
    d = d[d.argsort()[select]]
    if d < d_min:
      return False, d

  return True, -1

####################################################################################################

def check_structure(lattice, element, position, element_types, d_mins, path, timestep):
  '''
  '''
  N = len(element_types)
  for i in range(N):
    for j in range(i,N):
      accepted, d = check_nearest_neighbours(lattice, position[element==element_types[i]], position[element==element_types[j]], i==j, d_mins[element_types[i]][j-i])
      if not accepted:
        print('Too small interatomic distance in {0}_s{1}: {2}-{3}: {4} Ang'.format(path, timestep, element_types[i], element_types[j], d))
        return False

  return True

####################################################################################################

def read_structure(data, path, timestep, element_types, d_mins, names, lattices, elements, positions, charges):
  '''
  '''
  lattice, element, position, charge = get_structure(data)
  accepted = check_structure(lattice, element, position, element_types, d_mins, path, timestep)
  if accepted:
    names.append(path+'_s'+str(timestep))
    lattices.append(lattice)
    elements.append(element)
    positions.append(position)
    charges.append(charge)

  return accepted

####################################################################################################

def read_structures(path, extrapolation_data, selected_timestep, n_small, small, tolerance_index, extrapolation_statistic, element2index, element_types, d_mins, min_timestep_separation_extrapolation, names, lattices, elements, positions, charges):
  '''
  '''
  with open('mode1/'+path+'/structure.lammpstrj',encoding="utf-8") as f:
    data = [line.strip() for line in f.readlines()]

  n_interpolation_checks = len(selected_timestep)-2
  if n_interpolation_checks>0:
    for i in range(n_interpolation_checks):
      if selected_timestep[i]>=0:
        if selected_timestep[i]!=int(data[3]):
          start = data.index(str(selected_timestep[i]))-1
        else:
          start = 1
          while selected_timestep[i]!=int(data[start]):
            start += extrapolation_data[3]
          start -= 1
        end = start+extrapolation_data[3]
        accepted = read_structure(data[start:end], path, selected_timestep[i], element_types, d_mins, names, lattices, elements, positions, charges)
        if accepted:
          statistics.append([])
        data = data[end:]

  if selected_timestep[-2]>=min_timestep_separation_extrapolation:
    if selected_timestep[-2]!=int(data[3]):
      start = data.index(str(selected_timestep[-2]))-1
    else:
      start = 1
      while selected_timestep[-2]!=int(data[start]):
        start += extrapolation_data[3]
      start -= 1
    end = start+extrapolation_data[3]
    accepted = read_structure(data[start:end], path, selected_timestep[-2], element_types, d_mins, names, lattices, elements, positions, charges)
    if accepted:
      extrapolation_statistic[small] = np.array(extrapolation_statistic[small])
      extrapolation_statistic[small][:,1] = np.array([element2index[elements[-1][extrapolation_statistic[small][i,1].astype(int)-1]] for i in range(len(extrapolation_statistic[small]))])
      extrapolation_statistic[small] = extrapolation_statistic[small][extrapolation_statistic[small][:,0].argsort()]
      extrapolation_statistic[small] = extrapolation_statistic[small][extrapolation_statistic[small][:,2].argsort(kind='mergesort')]
      extrapolation_statistic[small] = extrapolation_statistic[small][extrapolation_statistic[small][:,1].argsort(kind='mergesort')]
      statistics.append(['small', str(list(np.array(element_types)[extrapolation_statistic[small][:,1].astype(int)])).strip('[]').replace("'",''), str(list(extrapolation_statistic[small][:,2].astype(int)+1)).strip('[]'), str([round(j, 5) for j in extrapolation_statistic[small][:,0]]).strip('[]')])

  accepted = False
  while not accepted and tolerance_index>=0:
    if selected_timestep[-1][tolerance_index]>=min_timestep_separation_extrapolation:
      if selected_timestep[-1][tolerance_index]!=int(data[3]):
        start = data.index(str(selected_timestep[-1][tolerance_index]))-1
      else:
        start = 1
        while selected_timestep[-1][tolerance_index]!=int(data[start]):
          start += extrapolation_data[3]
        start -= 1
      end = start+extrapolation_data[3]
      accepted = read_structure(data[start:end], path, selected_timestep[-1][tolerance_index], element_types, d_mins, names, lattices, elements, positions, charges)
    else:
      tolerance_index = -1
    if not accepted:
      tolerance_index -= 1
      if selected_timestep[-1][tolerance_index]-min_timestep_separation_extrapolation<selected_timestep[-2]:
        tolerance_index = -1
  if accepted:
    extrapolation_statistic[tolerance_index+n_small] = np.array(extrapolation_statistic[tolerance_index+n_small])
    extrapolation_statistic[tolerance_index+n_small][:,1] = np.array([element2index[elements[-1][extrapolation_statistic[tolerance_index+n_small][i,1].astype(int)-1]] for i in range(len(extrapolation_statistic[tolerance_index+n_small]))])
    extrapolation_statistic[tolerance_index+n_small] = extrapolation_statistic[tolerance_index+n_small][extrapolation_statistic[tolerance_index+n_small][:,0].argsort()]
    extrapolation_statistic[tolerance_index+n_small] = extrapolation_statistic[tolerance_index+n_small][extrapolation_statistic[tolerance_index+n_small][:,2].argsort(kind='mergesort')]
    extrapolation_statistic[tolerance_index+n_small] = extrapolation_statistic[tolerance_index+n_small][extrapolation_statistic[tolerance_index+n_small][:,1].argsort(kind='mergesort')]
    statistics.append(['large', str(list(np.array(element_types)[extrapolation_statistic[tolerance_index+n_small][:,1].astype(int)])).strip('[]').replace("'",''), str(list(extrapolation_statistic[tolerance_index+n_small][:,2].astype(int)+1)).strip('[]'), str([round(j, 5) for j in extrapolation_statistic[tolerance_index+n_small][:,0]]).strip('[]')])

  return tolerance_index

####################################################################################################

def write_data_new(names, lattices, elements, positions, charges, statistics, periodic):
  '''
  '''
  with open('input.data-new', 'w',encoding="utf-8") as f:
    for i in range(len(names)):
      f.write('begin\ncomment file {0}\n'.format(names[i]))
      if list(statistics[i]):
        f.write('comment statistics {0}\ncomment statistics {1}\ncomment statistics {2}\ncomment statistics {3}\n'.format(statistics[i][0], statistics[i][1], statistics[i][2], statistics[i][3]))
      if periodic:
        f.write('lattice {0:>9.5f} {1:>9.5f} {2:>9.5f}\nlattice {3:>9.5f} {4:>9.5f} {5:>9.5f}\nlattice {6:>9.5f} {7:>9.5f} {8:>9.5f}\n'.format(round(lattices[i][0][0]/Bohr2Ang,5), round(lattices[i][0][1]/Bohr2Ang,5), round(lattices[i][0][2]/Bohr2Ang,5), round(lattices[i][1][0]/Bohr2Ang,5), round(lattices[i][1][1]/Bohr2Ang,5), round(lattices[i][1][2]/Bohr2Ang,5), round(lattices[i][2][0]/Bohr2Ang,5), round(lattices[i][2][1]/Bohr2Ang,5), round(lattices[i][2][2]/Bohr2Ang,5)))
      for j in range(len(elements[i])):
        f.write('atom {0:>9.5f} {1:>9.5f} {2:>9.5f} {3:2} {4:>6.3f} 0.0 0.0 0.0 0.0\n'.format(round(positions[i][j][0]/Bohr2Ang,5), round(positions[i][j][1]/Bohr2Ang,5), round(positions[i][j][2]/Bohr2Ang,5), elements[i][j], charges[i][j]))
      f.write('energy 0.0\ncharge 0.0\nend\n')

####################################################################################################

def print_reliability(extrapolation_timesteps, smalls, tolerance_indices, paths):
  '''
  '''
  with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=RuntimeWarning)
    small_extrapolation_timesteps = np.diagonal(extrapolation_timesteps.T[smalls])
    extrapolation_timesteps_reduced = small_extrapolation_timesteps[small_extrapolation_timesteps!=-1]
    paths_reduced = paths[small_extrapolation_timesteps!=-1]
    median_small = np.median(extrapolation_timesteps_reduced)
    if np.isfinite(median_small):
      median_small = int(round(median_small,0))
      median_small_1 = np.median(extrapolation_timesteps_reduced[np.flatnonzero(np.core.defchararray.find(paths_reduced,'_hdnnp1_')!=-1)])
      if np.isfinite(median_small_1):
        median_small_1 = int(round(median_small_1,0))
      median_small_2 = np.median(extrapolation_timesteps_reduced[np.flatnonzero(np.core.defchararray.find(paths_reduced,'_hdnnp2_')!=-1)])
      if np.isfinite(median_small_2):
        median_small_2 = int(round(median_small_2,0))
      print('The median number of time steps to a small extrapolation is {0} (HDNNP_1: {1}, HDNNP_2: {2}).'.format(median_small, median_small_1, median_small_2))
    extra_steps = np.diagonal(extrapolation_timesteps.T[tolerance_indices])[tolerance_indices>=0]-small_extrapolation_timesteps[tolerance_indices>=0]
    if not np.isscalar(extra_steps):
      paths_reduced = paths[tolerance_indices>=0]
      median_extra_steps = np.median(extra_steps)
      if np.isfinite(median_extra_steps):
        median_extra_steps = int(round(median_extra_steps,0))
        median_extra_steps_1 = np.median(extra_steps[np.flatnonzero(np.core.defchararray.find(paths_reduced,'_hdnnp1_')!=-1)])
        if np.isfinite(median_extra_steps_1):
          median_extra_steps_1 = int(round(median_extra_steps_1,0))
        median_extra_steps_2 = np.median(extra_steps[np.flatnonzero(np.core.defchararray.find(paths_reduced,'_hdnnp2_')!=-1)])
        if np.isfinite(median_extra_steps_2):
          median_extra_steps_2 = int(round(median_extra_steps_2,0))
        print('The median number of time steps between the first and second selected extrapolated structure is {0} (HDNNP_1: {1}, HDNNP_2: {2}).'.format(median_extra_steps, median_extra_steps_1, median_extra_steps_2))

####################################################################################################

####################################################################################################

def read_data():
  '''
  '''
  names = []
  lattices = []
  elements = []
  positions = []
  charges = []
  statistics = []
  with open('input.data-new',encoding="utf-8") as f:
    for line in f.readlines():
      if line.startswith('atom'):
        line = line.split()
        elements[-1].append(line[4])
        positions[-1].append(line[1:4])
        charges[-1].append(line[5])
      elif line.startswith('lattice'):
        lattices[-1].append(line.strip().split()[1:4])
      elif line.startswith('comment file'):
        names.append(line.strip().split()[2])
      elif line.startswith('comment statistics'):
        statistics[-1].append(line.strip()[19:])
      elif line.startswith('begin'):
        lattices.append([])
        elements.append([])
        positions.append([])
        charges.append([])
        statistics.append([])
      elif line.startswith('end'):
        lattices[-1] = np.array(lattices[-1]).astype(float)*Bohr2Ang
        elements[-1] = np.array(elements[-1])
        positions[-1] = np.array(positions[-1]).astype(float)*Bohr2Ang
        charges[-1] = np.array(charges[-1]).astype(float)
  names = np.array(names)
  lattices = np.array(lattices)
  elements = np.array(elements)
  positions = np.array(positions)
  charges = np.array(charges)
  statistics = np.array(statistics)

  return names, lattices, elements, positions, charges, statistics

####################################################################################################

def print_performance(n_calculations):
  '''
  '''
  time = []
  for input_name in ['mode2/HDNNP_1/mode_1.out', 'mode2/HDNNP_1/mode_2.out', 'mode2/HDNNP_2/mode_1.out', 'mode2/HDNNP_2/mode_2.out']:
    with open(input_name,encoding="utf-8") as f:
      time.append([line.strip().split() for line in f.readlines() if line.startswith(' Total runtime (s)  :')][0])
  time = np.array(time)[:,4].astype(float)
  time = [time[0]+time[1], time[2]+time[3]]
  unit = ['s', 's']
  for i in range(2):
    if time[i]>=60.0:
      time[i] /= 60.0
      unit[i] = 'min'
      if time[i]>=60.0:
        time[i] /= 60.0
        unit[i] = 'h'
  print('\nTime to calculate {0} structures using RuNNer: HDNNP_1: {1} {2}, HDNNP_2: {3} {4}.\n'.format(n_calculations, round(time[0],2), unit[0], round(time[1],2), unit[1]))

####################################################################################################

def read_energies(input_name):
  '''
  '''
  with open(input_name,encoding="utf-8") as f:
    line = f.readline().strip()
    if line.startswith('point'):
      energies = np.array([line.strip().split()[2] for line in f.readlines()]).astype(float)
      energies = np.dstack((np.arange(len(energies)),energies))[0]
    elif line.startswith('Conf.'):
      energies = np.array([np.array(line.strip().split())[[1,3]] for line in f.readlines()]).astype(float)
      energies = energies[:,1]/energies[:,0]
      energies = np.dstack((np.arange(len(energies)),energies))[0]
    else:
      print('ERROR: Unknown RuNNer format')
      exit()

  return energies

####################################################################################################

def read_forces(input_name):
  '''
  '''
  with open(input_name,encoding="utf-8") as f:
    line = f.readline().strip()
    if line.startswith('point'):
      forces = np.array([np.array(line.strip().split())[[0,4]] for line in f if line.strip()]).astype(float)
      forces[:,0] -= 1
    elif line.startswith('Conf.'):
      forces = np.array([np.array(line.strip().split())[[0,5,6,7]] for line in f if line.strip()]).astype(float)
      forces = np.concatenate((forces[:,[0,1]], forces[:,[0,2]], forces[:,[0,3]]))
      forces[:,0] -= 1
    else:
      print('ERROR: Unknown RuNNer format')
      exit()

  return forces

####################################################################################################

def reduce_selection(selection, max_interpolated_structures_per_simulation, structure_name_index, timestep_separation_interpolation_checks, steps, indices):
  '''
  '''
  steps = np.array(steps)
  steps_difference = steps[1:]-steps[:-1]
  min_separation = steps_difference.min()
  if min_separation<timestep_separation_interpolation_checks[structure_name_index]:
    selection = selection[selection!=indices[1]]
  else:
    n_steps = len(steps)
    min_timestep_separation_interpolation_checks = (n_steps//max_interpolated_structures_per_simulation+1)*timestep_separation_interpolation_checks[structure_name_index]
    j = 1
    while j<n_steps-1:
      if steps_difference[j]<=min_timestep_separation_interpolation_checks:
        selection = selection[selection!=indices[j]]
        j += 1
      j += 1

  return selection

####################################################################################################

def improve_selection(selection, statistics, names, all_extrapolated_structures, max_extrapolated_structures, max_interpolated_structures_per_simulation, exceptions, structure_name_indices, timestep_separation_interpolation_checks):
  '''
  '''
  current_name = None
  steps = []
  indices = []
  for i in range(len(names)):
    if list(statistics[i]):
      if all_extrapolated_structures[structure_name_indices[i]]:
        selection = np.append(selection, i)
    elif i in selection:
      name = names[i].split('_s')
      if current_name==name[0]:
        steps.append(int(name[1]))
        indices.append(i)
      else:
        if len(steps)>2:
          selection = reduce_selection(selection, max_interpolated_structures_per_simulation[structure_name_indices[indices[0]]], structure_name_indices[indices[0]], timestep_separation_interpolation_checks, steps, indices)
        current_name = name[0]
        steps = [int(name[1])]
        indices = [i]
  if len(steps)>2:
    selection = reduce_selection(selection, max_interpolated_structures_per_simulation[structure_name_indices[indices[0]]], structure_name_indices[indices[0]], timestep_separation_interpolation_checks, steps, indices)
  selection = np.unique(selection)

  if any(max_extrapolated_structures) or any(exceptions):
    statistics_reduced = statistics[selection]
    structure_name_indices_reduced = structure_name_indices[selection]
    structure_name_indices_reduced = np.array([structure_name_indices_reduced[i] for i in range(len(structure_name_indices_reduced)) if list(statistics_reduced[i]) and statistics_reduced[i][0]=='small']).astype(str)
    if list(structure_name_indices_reduced):
      statistics_reduced = np.array([i for i in statistics_reduced if list(i) and i[0]=='small'])
      statistics_reduced = np.core.defchararray.add(np.core.defchararray.add(np.core.defchararray.add(np.core.defchararray.add(structure_name_indices_reduced, ';'), statistics_reduced[:,1]), ';'), statistics_reduced[:,2])
      statistics_unique = np.unique(statistics_reduced)
      counts = {}
      for i in statistics_unique:
        counts[i] = 0
      for i in statistics_reduced:
        counts[i] += 1
      exception_list = {}

      if any(max_extrapolated_structures):
        for i in statistics_unique:
          structure_index = int(i.split(';')[0])
          if max_extrapolated_structures[structure_index]!=0:
            if counts[i]>max_extrapolated_structures[structure_index]:
              exception_list[i] = np.concatenate((np.ones(max_extrapolated_structures[structure_index], dtype=int), np.zeros(counts[i]-max_extrapolated_structures[structure_index], dtype=int)))
              np.random.shuffle(exception_list[i])
              print('The extrapolation [\'{0}\', \'{1}\'] occurred {2} times.'.format(i.split(';')[1], i.split(';')[2], counts[i]))

      if any(exceptions):
        exceptions_unique = []
        for i in range(len(exceptions)):
          if exceptions[i]!=None:
            for j in range(len(exceptions[i])):
              exceptions_unique.append([str(i)+';'+exceptions[i][j][0]+';'+exceptions[i][j][1], exceptions[i][j][2]])
        counts_keys = counts.keys()
        for i in exceptions_unique:
          if i[0] in counts_keys:
            keep = int(round(i[1]*counts[i[0]]))
            exception_list[i[0]] = np.concatenate((np.ones(keep, dtype=int), np.zeros(counts[i[0]]-keep, dtype=int)))

      exception_list_keys = exception_list.keys()
      if list(exception_list_keys):
        structure_name_indices_reduced = structure_name_indices[selection]
        statistics_reduced = np.array(list(statistics[selection]))
        for i in range(len(selection)):
          if list(statistics_reduced[i]) and statistics_reduced[i][0]=='small':
            key = str(structure_name_indices_reduced[i])+';'+statistics_reduced[i][1]+';'+statistics_reduced[i][2]
            if key in exception_list_keys:
              if exception_list[key][-1]==0:
                selection[i] = -1
              exception_list[key] = np.delete(exception_list[key], -1, 0)
        selection = np.unique(selection)
        if selection[0] == -1:
          selection = selection[1:]

  return selection

####################################################################################################

def write_data_add(names, lattices, elements, positions, charges, periodic):
  '''
  '''
  with open('input.data-add', 'w',encoding="utf-8") as f:
    for i in range(len(names)):
      f.write('begin\ncomment file {0}\n'.format(names[i]))
      if periodic:
        f.write('lattice {0:>9.5f} {1:>9.5f} {2:>9.5f}\nlattice {3:>9.5f} {4:>9.5f} {5:>9.5f}\nlattice {6:>9.5f} {7:>9.5f} {8:>9.5f}\n'.format(round(lattices[i][0][0]/Bohr2Ang,5), round(lattices[i][0][1]/Bohr2Ang,5), round(lattices[i][0][2]/Bohr2Ang,5), round(lattices[i][1][0]/Bohr2Ang,5), round(lattices[i][1][1]/Bohr2Ang,5), round(lattices[i][1][2]/Bohr2Ang,5), round(lattices[i][2][0]/Bohr2Ang,5), round(lattices[i][2][1]/Bohr2Ang,5), round(lattices[i][2][2]/Bohr2Ang,5)))
      for j in range(len(elements[i])):
        f.write('atom {0:>9.5f} {1:>9.5f} {2:>9.5f} {3:2} {4:>6.3f} 0.0 0.0 0.0 0.0\n'.format(round(positions[i][j][0]/Bohr2Ang,5), round(positions[i][j][1]/Bohr2Ang,5), round(positions[i][j][2]/Bohr2Ang,5), elements[i][j], charges[i][j]))
      f.write('energy 0.0\ncharge 0.0\nend\n')

####################################################################################################

def print_statistics(selection, statistics, names, structure_names):
  '''
  '''
  if structure_names!=None and len(structure_names)>1:
    for structure_name in structure_names:
      print('Structure: {0}'.format(structure_name))
      n_extrapolations = int(np.array([1 for name in names[selection] if name.split('_')[0]==structure_name]).sum())
      print('{0} missing structures were identified.'.format(n_extrapolations))
      statistics_reduced = np.array([statistics[selection][i][0] for i in range(len(statistics[selection])) if names[selection][i].split('_')[0]==structure_name and list(statistics[selection][i])])
      if len(statistics_reduced)>0:
        n_small_extrapolations = int(np.array([1 for i in statistics_reduced if i=='small']).sum())
        n_large_extrapolations = int(np.array([1 for i in statistics_reduced if i=='large']).sum())
        print('{0} missing structures originate from small extrapolations.\n{1} missing structures originate from large extrapolations.'.format(n_small_extrapolations, n_large_extrapolations))
  else:
    print('{0} missing structures were identified.'.format(len(selection)))
    statistics_reduced = np.array([i[0] for i in statistics[selection] if list(i)])
    if len(statistics_reduced)>0:
      n_small_extrapolations = int(np.array([1 for i in statistics_reduced if i=='small']).sum())
      n_large_extrapolations = int(np.array([1 for i in statistics_reduced if i=='large']).sum())
      print('{0} missing structures originate from small extrapolations.\n{1} missing structures originate from large extrapolations.'.format(n_small_extrapolations, n_large_extrapolations))
  statistics = np.array([i for i in statistics[selection] if list(i)])
  if list(statistics):
    analyse_extrapolation_statistics(statistics)

####################################################################################################

def analyse_extrapolation_statistics(statistics):
  '''
  '''
  elements = []
  for line in statistics[:,1]:
    if ', ' in line:
      elements.extend(line.split(', '))
    else:
      elements.append(line)
  elements = np.array(elements)
  symmetry_functions = []
  for line in statistics[:,2]:
    if ', ' in line:
      symmetry_functions.extend(line.split(', '))
    else:
      symmetry_functions.append(line)
  symmetry_functions = np.array(symmetry_functions).astype(int)
  values = []
  for line in statistics[:,3]:
    if ', ' in line:
      values.extend(line.split(', '))
    else:
      values.append(line)
  values = np.array(values).astype(float)
  element_list = np.unique(elements)
  for e in element_list:
    symfunc = symmetry_functions[elements==e]
    symfunc_list = np.unique(symfunc)
    with open('extrapolation_statistics_'+e+'.dat', 'w',encoding="utf-8") as f:
      for s in symfunc_list:
        val = values[elements==e][symfunc==s]
        for v in val:
          f.write('{0} {1}\n'.format(s, v))

####################################################################################################

####################################################################################################

d_mins, min_timestep_separation_extrapolation, timestep_separation_interpolation_checks, min_timestep_separation_interpolation, delta_E, delta_F, all_extrapolated_structures, max_extrapolated_structures, max_interpolated_structures_per_simulation, exceptions = check_input(element_types, masses, structure_names, d_mins, min_timestep_separation_extrapolation, timestep_separation_interpolation_checks, min_timestep_separation_interpolation, delta_E, delta_F, all_extrapolated_structures, max_extrapolated_structures, max_interpolated_structures_per_simulation, exceptions, tolerances, initial_tolerance)

np.random.seed(seed)

restarted = False
if os.path.isfile('input.data-new'):
  print("The file input.data-new exists already and the data is reused. If the generation of the input.data-new file shall be redone, please remove the file.")
  restarted = True
  names, lattices, elements, positions, charges, statistics = read_data()
else:
  names = []
  lattices = []
  elements = []
  positions = []
  charges = []
  statistics = []
  if structure_names==None:
    n_structure_names = 1
  else:
    n_structure_names = len(structure_names)
  for i in range(n_structure_names):
    if structure_names==None:
      paths = get_paths('')
    else:
      print('Structure: {0}'.format(structure_names[i]))
      paths = get_paths(structure_names[i])
    extrapolation_data = []
    extrapolation_timesteps = []
    extrapolation_values = []
    extrapolation_statistics = []
    extrapolation_format = read_log_format(paths[0])
    for path in paths:
      extrapolation_data.append(read_extrapolation(path))
      extrapolation_timestep, extrapolation_value, extrapolation_statistic = read_log(path, extrapolation_data[-1], tolerances, extrapolation_format)
      extrapolation_timesteps.append(extrapolation_timestep)
      extrapolation_values.append(extrapolation_value)
      extrapolation_statistics.append(extrapolation_statistic)
    extrapolation_timesteps = np.array(extrapolation_timesteps).astype(int)
    extrapolation_values = np.array(extrapolation_values)
    selected_timesteps, tolerance_indices, smalls, n_small = get_timesteps(extrapolation_timesteps, extrapolation_values, extrapolation_data, min_timestep_separation_extrapolation[i], timestep_separation_interpolation_checks[i], min_timestep_separation_interpolation[i], tolerances, initial_tolerance)
    element2index = {}
    for j in range(len(element_types)):
      element2index[element_types[j]] = j
    for j in range(len(paths)):
      tolerance_indices[j] = read_structures(paths[j], extrapolation_data[j], selected_timesteps[j], n_small, smalls[j], tolerance_indices[j], extrapolation_statistics[j], element2index, element_types, d_mins[i], min_timestep_separation_extrapolation[i], names, lattices, elements, positions, charges)
    print_reliability(extrapolation_timesteps, smalls, tolerance_indices, paths)
  names = np.array(names)
  lattices = np.array(lattices)
  elements = np.array(elements)
  positions = np.array(positions)
  charges = np.array(charges)
  statistics = np.array(statistics)
  write_data_new(names, lattices, elements, positions, charges, statistics, periodic)

####################################################################################################

if os.path.isdir('mode2') and restarted:
  print("The directory mode2 exists already and the data is reused. If the RuNNer calculation shall be redone, please remove the directory.")
else:
  if os.path.isdir('mode2'):
    subprocess.Popen('rm -r mode2', shell=True).wait()
  command = "bash lib/RuNNerActiveLearn_2a.sh %s %s" % (str(delta_E[0]), str(delta_F[0]))
  subprocess.Popen(command.split(), shell=False).wait()

#print_performance(len(names))
#if structure_names==None:
#  structure_name_indices = np.zeros(len(names), dtype=int)
#else:
#  structure_name_indices = np.array([structure_names.index(name.split('_')[0]) for name in names])
#energies_1 = read_energies('mode2/HDNNP_1/trainpoints.000000.out')
#energies_2 = read_energies('mode2/HDNNP_2/trainpoints.000000.out')
#forces_1 = read_forces('mode2/HDNNP_1/trainforces.000000.out')
#forces_2 = read_forces('mode2/HDNNP_2/trainforces.000000.out')
#dE = np.array([delta_E[structure_name_index] for structure_name_index in structure_name_indices])
#dF = np.array([delta_F[structure_name_indices[i]] for i in range(len(structure_name_indices)) for j in range(3*len(positions[i]))])
#energies = energies_1[np.absolute(energies_2[:,1]-energies_1[:,1])>dE,0]
#forces = forces_1[np.absolute(forces_2[:,1]-forces_1[:,1])>dF,0]
#selection = np.unique(np.concatenate((energies, forces)).astype(int))
#selection = improve_selection(selection, statistics, names, all_extrapolated_structures, max_extrapolated_structures, max_interpolated_structures_per_simulation, exceptions, structure_name_indices, timestep_separation_interpolation_checks)
#write_data_add(names[selection], lattices[selection], elements[selection], positions[selection], charges[selection], periodic)
#print_statistics(selection, statistics, names, structure_names)
#
####################################################################################################

####################################################################################################
