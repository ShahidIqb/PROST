#!/usr/bin/env python
#from misc_yang import *
#sys.path.append('./source')

import time
import random
import argparse
import os, sys
import numpy as np
import numpy.matlib
import scipy.io as sp
from scipy import special
import load_bioinf_data as load_data


matmul = np.dot
def run_RNN(inputs, matrix_in, bias, mask):
  "Run an RNN layer for a given input, matrix of weights and a bias"
  # Init
  temp_x = None
  x = inputs
  output = np.zeros([inputs.shape[0], inputs.shape[1], matrix_in.shape[1]/4])
  h = np.matrix(np.zeros((inputs.shape[0], matrix_in.shape[1]/4)))
  c = np.matrix(np.zeros((inputs.shape[0], matrix_in.shape[1]/4)))
  #Loop over inputs
  for iterator in range(inputs.shape[1]):
      temp_x = x[:,iterator,:]
      temp_x = np.concatenate((temp_x, h), axis=1)
      inprod = matmul(temp_x, matrix_in)
      finalprod = inprod + bias

      # i = input_gate, j = new_input, f = forget_gate, o = output_gate
      i, j, f, o = np.split(finalprod, 4, axis=1)

      c = np.multiply(c, special.expit(f + 1.0)) + np.multiply(special.expit(i), np.tanh(j))
      h = np.multiply(np.tanh(c), special.expit(o))
      output[:,iterator,:] = h
  
  return (np.multiply(output, np.tile(mask[:,:,None], (1,1,matrix_in.shape[1]/4))), c)

def revseq(input, lengths):
    for i in range(input.shape[0]):
        input[i, 0:lengths[i], :] = input[i, (lengths[i]-1)::-1, :]
    return input

def rdelete(line, exts):
	""" return string by deleting 'ext' in the right """
	if isinstance(exts, str): exts = [exts]
	for ext in exts:
		if line.endswith(ext): return line[:-len(ext)]
	return line
#

def softmax(x):
    # """
    # Compute softmax values for each sets of scores in x.
    
    # Rows are scores for each class. 
    # Columns are predictions (samples).
    # """
    # scoreMatExp = np.exp(np.asarray(x))
    # return scoreMatExp / scoreMatExp.sum(0)
    a = np.exp(x)
    b = np.sum(np.exp(x), axis = 1)
    b = np.tile(b[:,None], (1,3))
    c = np.divide(a, b)
    return c
    #return np.exp(x) / np.sum(np.exp(x), axis=1)

def num_batches(num_vec, batch_size):
  incomplete_batch = 1 if np.mod(num_vec, batch_size) else 0
  return num_vec/batch_size+incomplete_batch

def bioinf_output_nonlinearity(output_type, pred):
  # this function applies the nonlinear activation functions for the different output types.
  
  output_type = output_type.upper()
  if output_type == 'SS':
    non_linear_output = softmax(pred)
  elif output_type == 'ASA' or output_type == 'HSEA' or output_type == 'HSEB' or output_type == 'CN' or output_type == 'CN13' or output_type == 'THETA' or output_type == 'TAU' or output_type == 'PHI' or output_type == 'PSI' or output_type == 'TT' or output_type == 'PP':
    non_linear_output = special.expit(pred)
  elif output_type == 'TTPP':
    non_linear_output = special.expit(pred)

  return non_linear_output
#
#
def run_all(list_dict_nn, list_pdb, opts):
	pssm_dir, hhm_dir, out_dir = opts.pssmdir, opts.hhmdir, opts.odir
	test_seq_names, test_feat0, test_aa = load_data.read_input(list_pdb, pssm_dir, hhm_dir)
	assert len(test_seq_names) == len(test_aa)
	test_lengths = [len(tmp) for tmp in test_feat0]

	cumlen1 = np.cumsum([0] + test_lengths)
	out12 = None

	for it1 in range(niter):
		dict_nn = list_dict_nn[it1]

		feat_mean, feat_var = dict_nn['feat']
		check_dat = dict_nn['wt']
		if it1 > 0:
			test_feat1 = [numpy.concatenate([a,b], axis=1) for a, b in zip(test_feat0, out12)]
		else:
			test_feat1 = test_feat0

		test_feat1 = load_data.do_mv_normalisation(test_feat1, None, input_mean=feat_mean, input_var=feat_var)[0]
		if opts.bdebug: print>>sys.stderr, 'printing ss', it1
		output_types = ['ss']
		out1 = run_batch1(check_dat[0], output_types, test_feat1)
		if opts.bdebug: print>>sys.stderr, 'printing the rest', it1
		output_types = ['asa', 'ttpp', 'hsea', 'hseb', 'cn']
		out2 = run_batch1(check_dat[1], output_types, test_feat1)

		out12 = np.concatenate([out1, out2], axis=1)
		out12 = [out12[a: cumlen1[k+1]] for k,a in enumerate(cumlen1[:-1])]

		if it1 != niter-1: continue

		assert len(out12) == len(test_seq_names) == len(test_aa)
		for pn, aa1, out1 in zip(test_seq_names, test_aa, out12):
			if it1 == niter-1:
				fout = os.path.join(out_dir, pn+opts.out_ext)
			else:
				fout = os.path.join(out_dir, pn+opts.out_ext+'_%d'%it1)
			load_data.write_spd33(fout, aa1, out1)
#		exit()
#	batch_seq_names = test_seq_names
#	batch_seq_lengths  = [len(x) for x in test_feat0]
#	misc_functions.save_predictions_to_file(out2, batch_seq_names, batch_seq_lengths, save_dir=outdir, file_ext='.i0r', header='%s' % ', '.join(map(str, output_types)))

def run_batch1(checkpoint_data, output_types, test_feat):

	true_label_ind, pred_label_ind, n_classes = load_data.get_outputs_list_stub(output_types)

	test_lengths = [len(tmp) for tmp in test_feat]

	# Network Parameters
	# n_input is the size of the features, i.e. 20 for PSSM
	n_input = len(test_feat[0][0])

	output_weights = checkpoint_data[0]
	output_bias = checkpoint_data[6]
	FC0_weights = checkpoint_data[13]
	FC0_bias = checkpoint_data[9]
	FC1_weights = checkpoint_data[1]
	FC1_bias = checkpoint_data[7]

	RNN1_FW_matrix = checkpoint_data[10]
	RNN1_BW_matrix = checkpoint_data[2]
	RNN1_FW_bias = checkpoint_data[3]
	RNN1_BW_bias = checkpoint_data[4]

	RNN2_FW_matrix = checkpoint_data[12]
	RNN2_BW_matrix = checkpoint_data[11]
	RNN2_FW_bias = checkpoint_data[8]
	RNN2_BW_bias = checkpoint_data[5]

	input_feat = test_feat
	batch_size = len(input_feat)+1	# NOW only one batch is permitted here
	seq_len = test_lengths

	for i in xrange(0, num_batches(len(input_feat), batch_size)):
		#print(i)
		#print "Doing batch ", i
		'''batch_ind = []
		for ib in range(int(i*batch_size), np.minimum(int((i+1)*batch_size), int(len(input_feat)))):
			batch_ind.append(ib)'''
		batch_ind = range(i*batch_size, np.minimum((i+1)*batch_size, len(input_feat)))
		#print(batch_ind)
		batch_seq_lengths = [ seq_len[ind] for ind in batch_ind ]
		batch_max_length = max(batch_seq_lengths)
		batch_feat = np.array( [ np.concatenate((np.array(tmp), np.zeros((batch_max_length - tmp.shape[0], len(input_feat[0][0]))))) for tmp in [ input_feat[ind] for ind in batch_ind ] ] )
		batch_seq_len_mask = np.array( [ np.concatenate((np.ones(tmp), np.zeros(batch_max_length - tmp))) for tmp in batch_seq_lengths ] )
#		batch_seq_names = [ test_seq_names[ind] for ind in batch_ind ]

		# Train/Test
		#Init and 0-pad
		#print(batch_seq_len_mask)
		x = batch_feat
		#print(x)
		# Run Forwards RNN1
		RNN1_fw, RNN1_state_fw = run_RNN(x, RNN1_FW_matrix, RNN1_FW_bias, batch_seq_len_mask)
		# Run Back RNN1 (reverse input and then reverse output)
		RNN1_bw, RNN1_state_bw = run_RNN(revseq(x, batch_seq_lengths), RNN1_BW_matrix, RNN1_BW_bias, batch_seq_len_mask)
		RNN1_bw = revseq(RNN1_bw, batch_seq_lengths)
		# Combine RNN1
		RNN1_out = np.concatenate((RNN1_fw, RNN1_bw), axis=2)
		# Run Forwards RNN2
		RNN2_fw, RNN2_state_fw = run_RNN(RNN1_out, RNN2_FW_matrix, RNN2_FW_bias, batch_seq_len_mask)
		# Run Backwards RNN2
		RNN2_bw, RNN2_state_bw = run_RNN(revseq(RNN1_out, batch_seq_lengths), RNN2_BW_matrix, RNN2_BW_bias, batch_seq_len_mask)
		RNN2_bw = revseq(RNN2_bw, batch_seq_lengths)
		# Combine RNN2
		RNN2_out = np.concatenate((RNN2_fw, RNN2_bw), axis=2)
		FC_in = RNN2_out[0, 0:batch_seq_lengths[0], :]
		for i in range(1,RNN2_out.shape[0]):
			FC_in = np.concatenate((FC_in, RNN2_out[i, 0:batch_seq_lengths[i], :]))

		# FC Layer 1
		FC0 = matmul(FC_in, FC0_weights)
		FC0 = np.add(FC0, FC0_bias)
		FC0_out = np.maximum(FC0, 0)
		# FC Layer 2
		FC1 = matmul(FC0_out, FC1_weights)
		FC1 = np.add(FC1, FC1_bias)
		FC1_out = np.maximum(FC1, 0)
		# Output Layer
		OUT = matmul(FC1_out, output_weights)
		OUT_out = np.add(OUT, output_bias)

		pred = OUT_out
		linear_output = pred
		
		output_index_pred = pred_label_ind

		temp_non_linear_output = []

		for ind, output_type in enumerate(output_types):
			temp_non_linear_output.append(bioinf_output_nonlinearity(output_type, linear_output[:, output_index_pred[ind][0]:output_index_pred[ind][1]]))
		
		non_linear_output = np.concatenate((temp_non_linear_output), axis=1)
		return non_linear_output
		#non_linear_output = softmax(pred)
		#We're done calculations. Save the results.
#		if not os.path.exists(outdir):
#			print "Making directory " + outdir
#			os.makedirs(outdir)

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('args', nargs='+', help='args')
	parser.add_argument('--nndir', dest="nndir",
						help="directory where the network and normalisation values are saved: default $xdir/../network_files/")
	parser.add_argument('--pssmdir', dest='pssmdir', default='.', help='pssm directory')
	parser.add_argument('--hhmdir', dest='hhmdir', default='.', help='hhm directory')
	parser.add_argument('--odir', dest='odir', default='.', help='out directory')
	parser.add_argument('--out_ext', dest='out_ext', default='.spd33', help='out file suffix')

	parser.add_argument("-f", dest='bforce', action='store_true', help='force mode')
	parser.add_argument("-v", dest='bdebug', action='store_true', help='debug mode')
	opts = parser.parse_args()
	if opts.nndir is None:
		opts.nndir = os.path.join(os.path.dirname(sys.argv[0]), '../network_files/')

	niter = 4
	list_dict_nn = None
	
	np.random.shuffle(opts.args)
	for i,x in enumerate(opts.args):
		bn = rdelete(os.path.basename(x), ['.seq', '.pssm', 'fa', '.hhm'])
		fout = os.path.join(opts.odir, bn+opts.out_ext)
		if not opts.bforce and os.path.isfile(fout): continue
		open(fout, 'w').close()

		if list_dict_nn is None:
			list_dict_nn = [np.load(opts.nndir + 'nndat_i%d.npz' % k, encoding='latin1') for k in range(niter)]

		run_all(list_dict_nn, [bn], opts)
