import os
#import h5py as h5
import numpy as np
import scipy.io as sio


###############################################################################
def load_casp_wrapper(casp_list_filename, feat_mean, feat_var):
  # this function will take a file as input where each line in the file has 3 fields
  # SEQ_NAME SEQ_PSSM_FILE_PATH SEQ_HMM_FILE_PATH
  # and there is one file for each of our casp sequences.
  # This is to be used for testing casp11 data, however can easilly be extended to any data.
  
  fp = open(casp_list_filename)
  lines = fp.readlines()
  fp.close()
  
  all_seq_names = []
  all_seq_data = []
  
  for line in lines:
    
    temp_line = line.split()
    
    all_seq_names.append(temp_line[0])
    seq_pssm = np.array(read_pssm(temp_line[1])[1])
    seq_hmm  = np.array(read_hmm(temp_line[2])[1])
    seq_data = np.concatenate([seq_pssm, seq_hmm], axis=1)
    
    all_seq_data.append(seq_data)
    
  # normalise the data
  all_seq_data_norm = do_mv_normalisation(all_seq_data, normalisation_mask=None, input_mean=feat_mean, input_var=feat_var)[0]
  
  return all_seq_names, all_seq_data_norm

###############################################################################  
def load_spd3_input_wrapper2(filename_list, feat_mean, feat_var, input_file_dir=None, input_file_ext='.spd3'):
	all_seq_names, all_seq_data = load_spd3_input_wrapper(filename_list, input_file_dir, input_file_ext)
	all_seq_data_norm = do_mv_normalisation(all_seq_data, None, input_mean=feat_mean, input_var=feat_var)[0]
	feature_length = len(all_seq_data[0][0])
	return all_seq_names, all_seq_data_norm, feature_length

def read_input(list_pdb, pssm_dir='.', hhm_dir='.', pssm_ext='.pssm'):
  # this function will take a file as input where each line in the file has 3 fields
  # SEQ_NAME SEQ_PSSM_FILE_PATH SEQ_HMM_FILE_PATH
  # we get the primary sequence from the pssm file.
  # we want to output pssm + hmm + phys7.
  
  # define the dictionary with the phys properties for each AA
  phys_dic = {'A': [-0.350, -0.680, -0.677, -0.171, -0.170, 0.900, -0.476],
              'C': [-0.140, -0.329, -0.359, 0.508, -0.114, -0.652, 0.476],
              'D': [-0.213, -0.417, -0.281, -0.767, -0.900, -0.155, -0.635],
              'E': [-0.230, -0.241, -0.058, -0.696, -0.868, 0.900, -0.582],
              'F': [ 0.363, 0.373, 0.412, 0.646, -0.272, 0.155, 0.318],
              'G': [-0.900, -0.900, -0.900, -0.342, -0.179, -0.900, -0.900],
              'H': [ 0.384, 0.110, 0.138, -0.271, 0.195, -0.031, -0.106],
              'I': [ 0.900, -0.066, -0.009, 0.652, -0.186, 0.155, 0.688],
              'K': [-0.088, 0.066, 0.163, -0.889, 0.727, 0.279, -0.265],
              'L': [ 0.213, -0.066, -0.009, 0.596, -0.186, 0.714, -0.053],
              'M': [ 0.110, 0.066, 0.087, 0.337, -0.262, 0.652, -0.001],
              'N': [-0.213, -0.329, -0.243, -0.674, -0.075, -0.403, -0.529],
              'P': [ 0.247, -0.900, -0.294, 0.055, -0.010, -0.900, 0.106],
              'Q': [-0.230, -0.110, -0.020, -0.464, -0.276, 0.528, -0.371],
              'R': [ 0.105, 0.373, 0.466, -0.900, 0.900, 0.528, -0.371],
              'S': [-0.337, -0.637, -0.544, -0.364, -0.265, -0.466, -0.212],
              'T': [ 0.402, -0.417, -0.321, -0.199, -0.288, -0.403, 0.212],
              'V': [ 0.677, -0.285, -0.232, 0.331, -0.191, -0.031, 0.900],
              'W': [ 0.479, 0.900, 0.900, 0.900, -0.209, 0.279, 0.529],
              'Y': [ 0.363, 0.417, 0.541, 0.188, -0.274, -0.155, 0.476],
              'X': [ 0.0771,-0.1536, -0.0620, -0.0762, -0.1451,  0.0497, -0.0398],
              'Z': [ 0.0771,-0.1536, -0.0620, -0.0762, -0.1451,  0.0497, -0.0398]}
  
  all_seq_names = []
  all_seq_data = []
  all_seq_aa = []
  
  for pdb1 in list_pdb:
    aa_pssm, pssm = read_pssm(os.path.join(pssm_dir, pdb1+pssm_ext))
    seq_phys = np.array( [ phys_dic[i] for i in aa_pssm] )
    seq_pssm = np.array(pssm)
    aa_hmm, seq_hmm = read_hmm( os.path.join(hhm_dir, pdb1+'.hhm'))
    if aa_hmm != aa_pssm:
      print>>stderr, pdb1, 'pssm/hhm seq not matching', len(aa_hmm), len(aa_pssm)
      if len(aa_hmm) == len(aa_pssm): print >>stderr, '\n'.join([aa_hmm, aa_pssm])
      continue
#
    seq_hmm  = np.array(seq_hmm)
    seq_data = np.concatenate([seq_pssm, seq_phys, seq_hmm], axis=1)
#    if input_file_dir is not None:
#      seq_file_input = load_spd3_file(input_file_dir + '/' + temp_line[0] + input_file_ext)
#      seq_data = np.concatenate([seq_pssm, seq_phys, seq_hmm, seq_file_input], axis=1)
#    else:
    
    all_seq_names.append(pdb1)
    all_seq_data.append(seq_data)
    all_seq_aa.append(aa_pssm)
    
  return all_seq_names, all_seq_data, all_seq_aa
 

###############################################################################
def read_pssm(pssm_file):
  # this function reads the pssm file given as input, and returns a LEN x 20 matrix (list) of pssm values.

  # index of 'ACDE..' in 'ARNDCQEGHILKMFPSTWYV'(blast order)
  idx_res = (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)
  
  # open the two files, read in their data and then close them
  fp = open(pssm_file, 'r')
  lines = fp.readlines()
  fp.close()

  # declare the empty dictionary with each of the entries
  aa = []
  pssm = []
  
  # iterate over the pssm file and get the needed information out
  for line in lines:
    split_line = line.replace('-', ' -').split()
    # valid lines should have 32 points of data.
    # any line starting with a # is ignored
    if (len(split_line) in (44,22)) and (split_line[0] not in ('#', 'Last')):
      aa_temp = split_line[1]
      aa.append(aa_temp)
      pssm_temp = [-float(i) for i in split_line[2:22]]
      pssm.append([pssm_temp[k] for k in idx_res])
  
  return aa, pssm


###############################################################################
def read_hmm(hhm_file):
  f = open(hhm_file)
  line=f.readline()
  while line[0]!='#':
      line=f.readline()
  f.readline()
  f.readline()
  f.readline()
  f.readline()
  seq = []
  extras = np.zeros([0,10])
  prob = np.zeros([0,20])
  line = f.readline()
  while line[0:2]!='//':
      lineinfo = line.split()
      seq.append(lineinfo[0])  
      probs_ = [2**(-float(lineinfo[i])/1000) if lineinfo[i]!='*' else 0. for i in range(2,22)]
      prob = np.concatenate((prob,np.matrix(probs_)),axis=0)
      
      line = f.readline()
      lineinfo = line.split()
      extras_ = [2**(-float(lineinfo[i])/1000) if lineinfo[i]!='*' else 0. for i in range(0,10)]
      extras = np.concatenate((extras,np.matrix(extras_)),axis=0)
      
      line = f.readline()
      assert len(line.strip())==0
      
      line = f.readline()
  #return (''.join(seq),prob,extras)
  return (seq,np.concatenate((prob,extras),axis=1))

###############################################################################
def read_mat(filename, field):
  return sio.loadmat(filename)[field] 
  

###############################################################################
def get_seq_lengths(data):
  # this function returns a list of the lengths of the sequences in data
  
  lengths = np.array([data['aa'][i].shape[0] for i in range(0,data['aa'].shape[0]) ], dtype='int64')

  return lengths
  

###############################################################################
def get_seq_names(data):
  # this function returns a list of sequence names for the input data.
  
  names = [ str(data[i]['name'][0][:-2]) for i in range(0, len(data)) ]
  
  return names

  
  
############################################################################### 
def get_min_max(data):
  # this functions takes the output of pad_list() and finds the min and max values that will be used to do the 0-1 normalisation.
  
  if type(data) is np.ndarray:
    input_min = np.min(np.min(data, axis=1), axis=0)
    input_max = np.max(np.max(data, axis=1), axis=0)
  elif type(data) is list:
    input_min = np.min([np.min(tmp, axis=0) for tmp in data], axis=0)
    input_max = np.max([np.max(tmp, axis=0) for tmp in data], axis=0)

  return input_min, input_max
  
  
############################################################################### 
def get_mean_variance(data):
  # this functions takes input data as a list and finds the mean and var values that will be used to do the 0 mean unit variance normalisation.
  
#  if type(data) is np.ndarray:
#    input_mean = np.mean(np.concatenate(data), axis=0)
#    input_var = np.var(np.concatenate(data), axis=0)
  if type(data) is list:
    input_mean = np.mean(np.concatenate(data), axis=0)
    input_var = np.var(np.concatenate(data), axis=0) 
#    input_mean = np.mean([np.min(tmp, axis=0) for tmp in data], axis=0)
#    input_var = np.var([np.max(tmp, axis=0) for tmp in data], axis=0)
  return input_mean, input_var
  
 
############################################################################### 
def do_zo_normalisation(data, input_min=None, input_max=None):
  # does normalisation between 0.05 and 0.95
  
  
  if normalisation_mask is None:
    normalisation_mask = np.ones(data[0].shape[1]) # THIS MAY NOT WORK FOR NPARRAY?
  
  if input_min is None:
    input_min, input_max = get_min_max(data)

  # do the masking
  input_min[normalisation_mask==0] = 0
  input_max[normalisation_mask==0] = 1
    
  # do the normalisation from 0 to 1
  if type(data) is np.array:
    normalised_data = (data - input_min) / (input_max - input_min)
  if type(data) is list:
    normalised_data = [(tmp - input_min) / (input_max-input_min) for tmp in data]
  
  # shift the nomalisation from 0.05 to 0.95 (to give room for new data, if need be)
  if type(data) is np.array:
    normalised_data = normalised_data * 0.9 + 0.05
  elif type(data) is list:
    temp = [tmp * 0.9 + 0.05 for tmp in normalised_data]
    normalised_data = temp
  
  return normalised_data, input_min, input_max
  
  
############################################################################### 
def do_mv_normalisation(data, normalisation_mask=None, input_mean=None, input_var=None):
  # does 0 mean unit variance normalisation
    
  if normalisation_mask is None:
    normalisation_mask = np.ones(data[0].shape[1]) # THIS MAY NOT WORK FOR NPARRAY?
  
  if input_mean is None:
    input_mean, input_var = get_mean_variance(data)

  # do the masking
  input_mean[normalisation_mask==0] = 0
  input_var[normalisation_mask==0] = 1
    
  # do the normalisation
  if type(data) is np.array:
    normalised_data = (data - input_mean) / np.sqrt(input_var)
  if type(data) is list:
    normalised_data = [(tmp - input_mean) / np.sqrt(input_var) for tmp in data]
  
   
  return normalised_data, input_mean, input_var


  
###############################################################################  
def pad_list(data_list, max_lengths):
  # this function will take the list output by get_inputs_list (or output) and will return a padded 3D numpy array
  # [num_seq, max_seq_len, input_feat_size] --> NOTE: the only size that is changed here is the second axis size.
  
  if max_lengths is None:
    # find the lengths + max length of the data
    lengths = [data_list[i].shape[0] for i in xrange(len(data_list))]
    max_lengths = max(lengths)
  
  # initialise the array to be all zeros
  padded_data = np.zeros((len(data_list), max_lengths, data_list[0][0].shape[0]), dtype=float)
  
  # put in our data
  for i, data in enumerate(data_list):
    padded_data[i,0:data.shape[0],:] = data
    
  return padded_data     



###############################################################################
def process_inputs_for_one_seq_one_type(data, input_type):
  # this function takes a single sequence data as input, and a single input_type
  # and will return the appropriate input data.
  #
  # 'PSSM' - returns the PSSM values for the sequence
  #
  # 'PHYS7' - will return 7 physio chemical properties of each residue
  #
  # 'ONEHOT' - will return a 1-hot representation of the sequence
  #
  # 'BLOSUM' - return the blosum values for the sequence
  
  input_type = input_type.upper()

  if input_type == 'PSSM':
    _input = data['pssm'].astype(float)
  elif input_type == 'PHYS7' or input_type == 'PHYS':
    _input = data['phys'].astype(float)
  elif input_type == 'HMM30':
    _input = 2**(-data['HHMprob'] / 1000.)
  elif input_type == 'ONEHOT':
    aa_dict = {'A':[1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'C':[0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'D':[0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'E':[0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'F':[0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'G':[0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'H':[0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'I':[0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'K':[0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'L':[0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'M':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
               'N':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.],
               'P':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.],
               'Q':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.],
               'R':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.],
               'S':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.],
               'T':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.],
               'V':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.],
               'W':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.],
               'Y':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.]}
    _input = np.array([aa_dict[i] for i in data['aa']])
  elif input_type == 'BLOSUM' or input_type == 'BLOSUM50':
    blosum_dict = {'A':[ 5.,-2.,-1.,-2.,-1.,-1.,-1., 0.,-2.,-1.,-2.,-1.,-1.,-3.,-1., 1., 0.,-3.,-2., 0.],
               'C':[-1.,-4.,-2.,-4.,13.,-3.,-3.,-3.,-3.,-2.,-2.,-3.,-2.,-2.,-4.,-1.,-1.,-5.,-3.,-1.],
               'D':[-2.,-2., 2., 8.,-4., 0., 2.,-1.,-1.,-4.,-4.,-1.,-4.,-5.,-1., 0.,-1.,-5.,-3.,-4.],
               'E':[-1., 0., 0., 2.,-3., 2., 6.,-3., 0.,-4.,-3., 1.,-2.,-3.,-1.,-1.,-1.,-3.,-2.,-3.],
               'F':[-3.,-3.,-4.,-5.,-2.,-4.,-3.,-4.,-1., 0., 1.,-4., 0., 8.,-4.,-3.,-2., 1., 4.,-1.],
               'G':[ 0.,-3., 0.,-1.,-3.,-2.,-3., 8.,-2.,-4.,-4.,-2.,-3.,-4.,-2., 0.,-2.,-3.,-3.,-4.],
               'H':[-2., 0., 1.,-1.,-3., 1., 0.,-2.,10.,-4.,-3., 0.,-1.,-1.,-2.,-1.,-2.,-3., 2.,-4.],
               'I':[-1.,-4.,-3.,-4.,-2.,-3.,-4.,-4.,-4., 5., 2.,-3., 2., 0.,-3.,-3.,-1.,-3.,-1., 4.],
               'K':[-1., 3., 0.,-1.,-3., 2., 1.,-2., 0.,-3.,-3., 6.,-2.,-4.,-1., 0.,-1.,-3.,-2.,-3.],
               'L':[-2.,-3.,-4.,-4.,-2.,-2.,-3.,-4.,-3., 2., 5.,-3., 3., 1.,-4.,-3.,-1.,-2.,-1., 1.],
               'M':[-1.,-2.,-2.,-4.,-2., 0.,-2.,-3.,-1., 2., 3.,-2., 7., 0.,-3.,-2.,-1.,-1., 0., 1.],
               'N':[-1.,-1., 7., 2.,-2., 0., 0., 0., 1.,-3 ,-4., 0.,-2.,-4.,-2., 1., 0.,-4.,-2.,-3.],
               'P':[-1.,-3.,-2.,-1.,-4.,-1.,-1.,-2.,-2.,-3.,-4.,-1.,-3.,-4.,10.,-1.,-1.,-4.,-3.,-3.],
               'Q':[-1., 1., 0., 0.,-3., 7., 2.,-2., 1.,-3.,-2., 2., 0.,-4.,-1., 0.,-1.,-1.,-1.,-3.],
               'R':[-2., 7.,-1.,-2.,-4., 1., 0.,-3., 0.,-4.,-3., 3.,-2.,-3.,-3.,-1.,-1.,-3.,-1.,-3.],
               'S':[ 1.,-1., 1., 0.,-1., 0.,-1., 0.,-1.,-3.,-3., 0.,-2.,-3.,-1., 5., 2.,-4.,-2.,-2.],
               'T':[ 0.,-1., 0.,-1.,-1.,-1.,-1.,-2.,-2.,-1.,-1.,-1.,-1.,-2.,-1., 2., 5.,-3.,-2., 0.],
               'V':[ 0.,-3.,-3.,-4.,-1.,-3.,-3.,-4.,-4., 4., 1.,-3., 1.,-1.,-3.,-2., 0.,-3.,-1., 5.],
               'W':[-3.,-3.,-4.,-5.,-5.,-1.,-3.,-3.,-3.,-3.,-2.,-3.,-1., 1.,-4.,-4.,-3.,15., 2.,-3.],
               'Y':[-2.,-1.,-2.,-3.,-3.,-1.,-2.,-3., 2.,-1.,-1.,-2., 0., 4.,-3.,-2.,-2., 2., 8.,-1.]}
    _input = np.array([blosum_dict[i] for i in data['aa']])
  else:    
    _input = None
    print ("ERROR GETTING INPUT DATA", input_type, "IS NOT VALID")

  return _input
  

###############################################################################  
def get_normalisation_mask_all_types(input_types):
  
  normalisation_mask = get_normalisation_mask_one_type(input_types[0])
  for input_type in input_types[1:]:
    normalisation_mask = np.concatenate([normalisation_mask, get_normalisation_mask_one_type(input_type)])
    
  
  return normalisation_mask
  
  
###############################################################################  
def get_normalisation_mask_one_type(input_type):

  input_type = input_type.upper()

  if input_type == 'PSSM':
    normalisation_mask = np.ones(20)
  elif input_type == 'PHYS7' or input_type == 'PHYS':
    normalisation_mask = np.ones(7)
  elif input_type == 'HMM30':
    normalisation_mask = np.ones(30)
  elif input_type == 'ONEHOT':
    normalisation_mask = np.zeros(20)
  elif input_type == 'BLOSUM' or input_type == 'BLOSUM50':
    normalisation_mask = np.ones(20)
  else:    
    normalisation_mask = None
    print ("ERROR GETTING INPUT DATA", input_type, "IS NOT VALID")
      
  return normalisation_mask

###############################################################################
def process_inputs_for_one_seq_all_types(data, input_types):
  # this function takes a single input sequence (data), and all of the output_types
  # and will return a numpy array of the output labels of shape [seq_len, sum(feat_len)]
  
  _input = process_inputs_for_one_seq_one_type(data, input_types[0])
  for _type in input_types[1:]:
    _temp_input = process_inputs_for_one_seq_one_type(data,_type)
    _input = np.concatenate((_input, _temp_input), axis=1)
    
    
  return _input


##############################################################################
def load_spd3_file(filename):
  # this function should load one file.
  
  read_data = np.loadtxt(filename)
  
  return read_data

##############################################################################
def load_spd3_files(data, input_file_dir=['./'], file_ext='.spd3'):
  # this function should load the files in the input_file_dir
  # data should be a list of the data loaded from the .mat files.
  # the function should return a list of the file contents, one element for each of the sequences in data.
  
  # get the names of the sequences
  names = get_seq_names(data)
  
  # load the data for each of the sequences
  read_data=[]
  for name in names:
    temp_read_data = []
    for input_dir in input_file_dir:
      str_name = input_dir+'/'+name+file_ext
      # I should really put some sort of test here to see if the file exists...
      temp_read_data.append(np.loadtxt(str_name))
    read_data.append(np.concatenate(temp_read_data, axis=1))
    
    
  return read_data
  

###############################################################################
def get_inputs_list(data, input_types):
  # this function takes all input sequences (data), and all of the output_types
  # and will return a numpy array of the output labels of shape [seq_len, sum(feat_len)]
  
  _input = [process_inputs_for_one_seq_all_types(_data, input_types) for _data in data]
  normalisation_mask = get_normalisation_mask_all_types(input_types)
  
  return _input, normalisation_mask
  
  
###############################################################################
def get_zo_normalised_inputs_list(data, input_types, input_min=None, input_max=None):
  # this function takes all input sequences (data), and all of the output_types
  # and will return a numpy array of the output labels of shape [seq_len, sum(feat_len)]
  
  input_list, normalisation_mask = get_inputs_list(data, input_types)
  
  # load data from text files.
  if input_file_dir is not None:
    read_data = load_spd3_files(data, input_file_dir)
    for ind, temp in enumerate(input_list):
      input_list[ind] = np.concatenate((input_list[ind], read_data[ind]), axis=1)    
    normalisation_mask = np.concatenate((normalisation_mask, np.ones(read_data[0].shape[1])))
    
  normalised_inputs, input_min, input_max = do_zo_normalisation(input_list, normalisation_mask, input_min, input_max)
  
  input_size = len(normalised_inputs[0][0])
  
  return normalised_inputs, input_min, input_max, input_size
  
  
###############################################################################
def get_mv_normalised_inputs_list(data, input_types, input_mean=None, input_var=None, input_file_dir=None, input_file_dir_ext='.spd3'):
  # this function takes all input sequences (data), and all of the output_types
  # and will return a numpy array of the output labels of shape [seq_len, sum(feat_len)]
      
  input_list, normalisation_mask = get_inputs_list(data, input_types)
  
  # load data from text files.
  if input_file_dir is not None:
    read_data = load_spd3_files(data, input_file_dir, file_ext=input_file_dir_ext)
    for ind, temp in enumerate(input_list):
      input_list[ind] = np.concatenate((input_list[ind], read_data[ind]), axis=1)    
    normalisation_mask = np.concatenate((normalisation_mask, np.ones(read_data[0].shape[1])))
    
  normalised_inputs, input_mean, input_var = do_mv_normalisation(input_list, normalisation_mask, input_mean, input_var)
  
  input_size = len(normalised_inputs[0][0])
  
  return normalised_inputs, input_mean, input_var, input_size

  
  
###############################################################################
def get_inputs(data, input_types, max_len=None):
  # this is the wrapper function that will simply return the numpy array of padded inputs
  # I want my inputs to be in the form of 3d np array [num_seq, max_seq_len, input_feat_size]
  # inputs: - data will be the output of our read_mat function
  #         - input_type will be a list of input types
  
  input_list = get_inputs_list(data, input_types)
  
  padded_inputs = pad_list(input_list, max_len)
  
  return padded_inputs
 
 

###############################################################################
def get_zo_normalised_inputs(data, input_types, max_len=None, input_min=None, input_max=None):
  
  input_list = get_inputs_list(data, input_types)
  
  padded_inputs = pad_list(input_list, max_len)
  
  normalised_inputs, input_min, input_max = do_zo_normalisation(padded_inputs, input_min, input_max)
  
  input_size = normalised_inputs.shape[2]
  
  return normalised_inputs, input_min, input_max, input_size
 
 

###############################################################################
def get_mv_normalised_inputs(data, input_types, max_len=None, input_mean=None, input_var=None):
  
  input_list = get_inputs_list(data, input_types)
  
  padded_inputs = pad_list(input_list, max_len)
  
  normalised_inputs, input_min, input_max = do_mv_normalisation(padded_inputs, input_mean, input_var)
  
  input_size = normalised_inputs.shape[2]
  
  return normalised_inputs, input_mean, input_var, input_size
  
  
  
###############################################################################  
def process_label_for_one_seq_one_type(data=None, output_type=None):
  # this function takes a single input sequence (data), and a single output_type 
  # and will return a numpy array of the output labels of shape [seq_len, feat_len].
  # this function will also find a mask for valid outputs. (1 for valid, 0 for invalid).
  # this mask can then be applied at the output (for example to the loss) so that the network
  # isn't trying to learn invalid outputs.
  #
  # 'SS' - will give a 1-hot array for the 3 state secondary structures 
  #        (will return all 0s for any X case) 
  #
  # 'ASA' - will return a 1D array of the ASA values (ASA/100)
  
  output_type = output_type.upper()
  
  if output_type == 'SS':
    if data is not None:
      ss_dict_one_hot = {"C":[1.,0.,0.], "H":[0.,1.,0.], "E":[0.,0.,1.], "X":[0.,0.,0.]}
      ss_dict = {"C":[0],"H":[1],"E":[2],"X":[-1]}
      ss_mask_dict = {"C":[1],"H":[1],"E":[1],"X":[0]}
      ss_mask_dict_one_hot = {"C":[1.,1.,1.], "H":[1.,1.,1.], "E":[1.,1.,1.], "X":[0.,0.,0.]}
      label = np.array([ ss_dict[i] for i in data['ss'] ])
      label_encoded = np.array([ ss_dict_one_hot[i] for i in data['ss'] ])
      mask = np.array([ ss_mask_dict[i] for i in data['ss'] ])
      mask_encoded = np.array([ ss_mask_dict_one_hot[i] for i in data['ss'] ])
    true_label_size = 1
    pred_label_size = 3
  elif output_type == 'ASA':
    if data is not None:
      label = np.array(data['asa']) / 100.0
      mask = np.array((data['asa']!=360).astype(float))
      # apply the mask here so that the 360 values aren't learned - shouldn't be needed
  #    label = label*mask
      label_encoded = label
      mask_encoded = mask
    true_label_size = 1
    pred_label_size = 1
  elif output_type == 'TTPP':  
    if data is not None:
      label = np.concatenate([data['t'], data['phi'], data['psi']], axis=1)
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.    
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      # apply the mask here so that the 360 values aren't learned - shouldn't be needed
  #    label = label*mask
      label_encoded = label
      mask_encoded = mask
    true_label_size = 8
    pred_label_size = 8
  elif output_type == 'THETA':  
    if data is not None:
      label = np.array([data['t'][:,0]]).T
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.    
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      label_encoded = label
      mask_encoded = mask
    true_label_size = 2
    pred_label_size = 2
  elif output_type == 'TAU':  
    if data is not None:
      label = np.array([data['t'][:,1]]).T
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.    
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      label_encoded = label
      mask_encoded = mask
    true_label_size = 2
    pred_label_size = 2
  elif output_type == 'TT':  
    if data is not None:
      label = data['t']
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      label_encoded = label
      mask_encoded = mask
    true_label_size = 4
    pred_label_size = 4
  elif output_type == 'PHI':  
    if data is not None:
      label = data['phi']
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.    
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      label_encoded = label
      mask_encoded = mask
    true_label_size = 2
    pred_label_size = 2
  elif output_type == 'PSI':  
    if data is not None:
      label = data['psi']
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.    
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      label_encoded = label
      mask_encoded = mask
    true_label_size = 2
    pred_label_size = 2
  elif output_type == 'PP':  
    if data is not None:
      label = np.concatenate([data['phi'], data['psi']], axis=1)
      mask = np.array(label!=360).astype(float)
      label = np.radians(label)
  #    label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)) * 0.99  # scale to -0.99 to 0.99 range. tanh activation.
      label = (np.concatenate([np.sin(label), np.cos(label)], axis=1)+1) / 2.  # scale to 0 to 1 range. sigmoid activation.
      mask = np.concatenate([mask, mask], axis=1)
      label_encoded = label
      mask_encoded = mask
    true_label_size = 4
    pred_label_size = 4
  elif output_type == 'HSEA':
    if data is not None:
      label = np.array(data['HSEa']) / [50., 65.]
      mask = (np.isnan(label).astype(float)==0).astype(float)
      label[np.isnan(label)] = 0 # set NaN values to be 0.
      label_encoded = label
      mask_encoded = mask
    true_label_size = 2
    pred_label_size = 2
  elif output_type == 'HSEB':
    if data is not None:
      label = np.array(data['HSEb']) / [50., 65.]
      mask = (np.isnan(label).astype(float)==0).astype(float)
      label[np.isnan(label)] = 0 # set NaN values to be 0.
      label_encoded = label
      mask_encoded = mask
    true_label_size = 2
    pred_label_size = 2
  elif output_type == 'CN' or output_type == 'CN13':
    if data is not None:
      label = np.array(data['CN13']) / 85.
      mask = (np.isnan(label).astype(float)==0).astype(float)
      label[np.isnan(label)] = 0 # set NaN values to be 0.
      label_encoded = label
      mask_encoded = mask
    true_label_size = 1
    pred_label_size = 1
  else:    
    label = None
    print ("ERROR GETTING OUTPUT LABELS ", output_type, " IS NOT VALID")

  if data is not None:
    return label, label_encoded, mask, mask_encoded, true_label_size, pred_label_size
  else:
    return None, None, None, None, true_label_size, pred_label_size
  



###############################################################################
def process_labels_for_one_seq_all_types(data=None, output_types=None):
  # this function takes a single input sequence (data), and all of the output_types
  # and will return a numpy array of the output labels of shape [seq_len, sum(feat_len)]
  
  
  labels, labels_encoded, mask, mask_encoded, true_label_size, pred_label_size = process_label_for_one_seq_one_type(data, output_types[0])
  list_of_true_label_sizes = [true_label_size]
  list_of_pred_label_sizes = [pred_label_size]
  
  for _type in output_types[1:]:
    t_labels, t_labels_encoded, t_mask, t_mask_encoded, t_true_label_size, t_pred_label_size = process_label_for_one_seq_one_type(data, _type)
    if data is not None:
      labels = np.concatenate((labels, t_labels), axis=1)
      labels_encoded = np.concatenate((labels_encoded, t_labels_encoded), axis=1)
      mask = np.concatenate((mask, t_mask), axis=1)
      mask_encoded = np.concatenate((mask_encoded, t_mask_encoded), axis=1)
    list_of_true_label_sizes = list_of_true_label_sizes + [t_true_label_size]   
    list_of_pred_label_sizes = list_of_pred_label_sizes + [t_pred_label_size]  
  
  if data is not None:  
    return labels, labels_encoded, mask, mask_encoded, list_of_true_label_sizes, list_of_pred_label_sizes
  else:
    return list_of_true_label_sizes, list_of_pred_label_sizes


###############################################################################
def get_outputs_list_stub(output_types):

  list_of_true_label_sizes, list_of_pred_label_sizes = process_labels_for_one_seq_all_types(output_types=output_types)
    
  list_of_true_label_sizes = [0] + list_of_true_label_sizes
  true_label_ind = [ [sum(list_of_true_label_sizes[0:i])] + [sum(list_of_true_label_sizes[0:i+1])] for i in range(1, len(list_of_true_label_sizes)) ]

  list_of_pred_label_sizes = [0] + list_of_pred_label_sizes  
  pred_label_ind = [ [sum(list_of_pred_label_sizes[0:i])] + [sum(list_of_pred_label_sizes[0:i+1])] for i in range(1, len(list_of_pred_label_sizes)) ]

  n_classes = sum(list_of_pred_label_sizes)  
  
  return true_label_ind, pred_label_ind, n_classes
  

###############################################################################
def get_outputs_list(data, output_types):
  # this function takes all input sequences (data), and all of the output_types
  # and will return a numpy array of the output labels of shape [seq_len, sum(feat_len)]
  
  _tuple = [process_labels_for_one_seq_all_types(_data, output_types) for _data in data]
  labels = [temp[0][:] for temp in _tuple]
  labels_encoded = [temp[1][:] for temp in _tuple]
  mask = [temp[2][:] for temp in _tuple]
  mask_encoded = [temp[3][:] for temp in _tuple]
  list_of_true_label_sizes = _tuple[0][4]
  list_of_pred_label_sizes = _tuple[0][5]
  n_classes = sum(list_of_pred_label_sizes)
    
  list_of_true_label_sizes = [0] + list_of_true_label_sizes
  true_label_ind = [ [sum(list_of_true_label_sizes[0:i])] + [sum(list_of_true_label_sizes[0:i+1])] for i in range(1, len(list_of_true_label_sizes)) ]

  list_of_pred_label_sizes = [0] + list_of_pred_label_sizes  
  pred_label_ind = [ [sum(list_of_pred_label_sizes[0:i])] + [sum(list_of_pred_label_sizes[0:i+1])] for i in range(1, len(list_of_pred_label_sizes)) ]
  
  return labels, labels_encoded, mask, mask_encoded, true_label_ind, pred_label_ind, n_classes
    

###############################################################################  
def get_outputs(data, output_types, max_len = None):
  # this is the wrapper function that will simply return the numpy array of padded outputs
  # I want my outputs to be in the form of 3d np array [num_seq, max_seq_len, output_label_size]
  # inputs: - data will be the output of our read_mat function
  #         - output_type will be a list of output types
  
  label_list, mask_list = get_outputs_list(data, output_types)
  
  padded_label = pad_list(label_list, max_len)
  padded_mask = pad_list(mask_list, max_len)
  
  output_size = padded_label.shape[2]
  
  return padded_label, padded_mask, output_size
  
#
def write_spd33(output_filename, seq1, raw_data):
  assert len(seq1) == len(raw_data)
  rnam1_std = "ACDEFGHIKLMNPQRSTVWY"
  ASA_std = (115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
			185, 160, 145, 180, 225, 115, 140, 155, 255, 230)
  dict_rnam1_ASA = dict(zip(rnam1_std, ASA_std))
  ASA0 = np.asarray([dict_rnam1_ASA.get(x, ASA_std[0]) for x in seq1])

  ss_order = ['C','H','E']
  ss_ind = np.argmax(raw_data[:,0:3], axis=1)
  pred_ss = np.array([ss_order[i] for i in ss_ind])

  pred_asa = raw_data[:,3] * ASA0

  raw_ttpp = raw_data[:,4:12] * 2 - 1;
  pred_theta = np.rad2deg(np.arctan2(raw_ttpp[:,0], raw_ttpp[:,4]))
  pred_tau = np.rad2deg(np.arctan2(raw_ttpp[:,1], raw_ttpp[:,5]))
  pred_phi = np.rad2deg(np.arctan2(raw_ttpp[:,2], raw_ttpp[:,6]))
  pred_psi = np.rad2deg(np.arctan2(raw_ttpp[:,3], raw_ttpp[:,7]))

  pred_hsea_up = raw_data[:,12] * 50.
  pred_hsea_down = raw_data[:,13] * 65.
  pred_hseb_up = raw_data[:,14] * 50.
  pred_hseb_down = raw_data[:,15] * 65.

  pred_cn = raw_data[:,16] * 85.

  readable_data = np.zeros(pred_ss.size,
          dtype=[
                 ('index', int),
                 ('pred_seq', 'S1'),
                 ('pred_ss', 'S1'),
                 ('pred_asa', float),
                 ('pred_phi', float),
                 ('pred_psi', float),
                 ('pred_theta', float),
                 ('pred_tau', float),
                 ('pred_hseau', float),
                 ('pred_hsead', float),
                 ('pred_pc', float),
                 ('pred_ph', float),
                 ('pred_pe', float) ])

  readable_data['index'] = np.arange(len(pred_ss)) + 1
  readable_data['pred_ss'] = pred_ss
  readable_data['pred_seq'] = np.array(list(seq1))
  readable_data['pred_asa'] = pred_asa
  readable_data['pred_phi'] = pred_phi
  readable_data['pred_psi'] = pred_psi
  readable_data['pred_theta'] = pred_theta
  readable_data['pred_tau'] = pred_tau
  readable_data['pred_hseau'] = pred_hsea_up
  readable_data['pred_hsead'] = pred_hsea_down
  readable_data['pred_pc'] = raw_data[:,0]
  readable_data['pred_ph'] = raw_data[:,1]
  readable_data['pred_pe'] = raw_data[:,2]
#  readable_data['pred_hsebu'] = pred_hseb_up
#  readable_data['pred_hsebd'] = pred_hseb_down
#  readable_data['pred_cn'] = pred_cn

  np.savetxt(output_filename, readable_data, fmt="%-3d %s %s %5.1f %6.1f %6.1f %6.1f %6.1f %4.1f %4.1f %5.3f %5.3f %5.3f", header='SEQ SS ASA Phi Psi Theta(i-1=>i+1) Tau(i-2=>i+2) HSE_alpha_up HSE_alpha_down P(C) P(H) P(E)')
