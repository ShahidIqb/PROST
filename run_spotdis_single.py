#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
#import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import tqdm,sys,glob,os
import numpy as np
#import cPickle as pickle
import _pickle as pickle
import pickle
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--gpu', default=-1, type=int, help='GPU number to use (default -1)')
parser.add_argument('--batch_size', default=50, type=int, help='Batch size to process inputs. Lower this value if running you''re out of memory.')
parser.add_argument('--output_dir', default='', type=str, help='Directory to place outputs in. If it''s unused then it places it where the input fasta file is.')
parser.add_argument('--quiet', default=0, type=int, help='Quiet (0=No, 1=Yes)')
parser.add_argument('vars', nargs='*')
args = parser.parse_args()


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)if args.gpu!=-1 else ""
os.environ["CUDA_VISIBLE_DEVICES"]="0"

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_dict = {}
for i in AA:
    one_hot = np.zeros([len(AA)])
    one_hot[AA.index(i)] = 1
    AA_dict[i] = one_hot

def read_fasta_file(fname):
    with open(fname,'r') as f:
        AA = ''.join(f.read().splitlines()[1:])
    return AA

def read_seq_file(fname):
    with open(fname,'r') as f:
        AA = f.read()
    return AA

def read_one_hot(fname):
    if fname[-3:] == '.seq':
        seq = read_seq_file(fname)
    else:
        seq = read_fasta_file(fname)
    return np.array([AA_dict.get(i,np.zeros([len(AA)])) for i in seq])

def sigmoid(x):
    return 1/(1+np.exp(-x))


'''with open('dat/norm_params.p','r') as f:
    normdic = pickle.load(f)'''
#normdic = pickle.load(open('dat/norm_params.p','r'))
with open('dat/norm_params.p', 'rb') as pickle_file:
    normdic = pickle.load(pickle_file, encoding='latin1')
normmu = normdic['mu']
normstd = normdic['std']

models = glob.glob('./dat/*meta')
all_outputs = []

if args.output_dir == '':
    base_name = [os.path.splitext(i)[0] for i in args.vars]
else:
    base_name = [args.output_dir + '/' + os.path.splitext(os.path.basename(i))[0] for i in args.vars]

print(args.vars)

tqdmfile = sys.stdout if args.quiet == 0 else os.devnull
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n/10%10!=1)*(n%10<4)*n%10::4])
for I,i in enumerate(models):
    model_id = i[:-5]
    if args.quiet == 0:        
        print('Going through %s model...'%(ordinal(I+1)))
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.allow_soft_placement = True
    config.log_device_placement = False
    saver = tf.train.import_meta_graph(i)
    with tf.Session() as sess:
        saver.restore(sess,model_id)
        output = tf.get_collection('output')[0]
        batch_size = args.batch_size
        num_batches = len(args.vars)/float(batch_size)
        tmp_outputs = []
        for j in tqdm.tqdm(np.arange(num_batches),file=tqdmfile):
            data = [(read_one_hot(args.vars[int(k)])-normmu)/normstd for k in np.arange(j*batch_size,min((j+1)*batch_size,len(args.vars)))]
            seq_lens = np.array([k.shape[0] for k in data])
            max_seq_len = np.max(seq_lens)
            mask = np.array([np.concatenate([np.ones([k]),np.zeros([max_seq_len-k])],0) for k in seq_lens])
            data = np.array([np.concatenate([k,np.zeros([max_seq_len-seq_lens[K],k.shape[1]])],0) for K,k in enumerate(data)])
            feed_dict = {'oneD_feats:0':data,'seq_lens:0':seq_lens,'ph_dropout:0':1.0,'mask_bool:0':mask,'train_bool:0':False,'ln_mask_bool:0':mask}
            outputs = sigmoid(sess.run(output,feed_dict=feed_dict))
            cum_seq_lens = np.concatenate([[0],np.cumsum(seq_lens)])
            tmp_outputs += [outputs[cum_seq_lens[k]:cum_seq_lens[k+1]] for k in range(len(cum_seq_lens)-1)]
        all_outputs.append({base_name[I]:tmp_outputs[I] for I,i in enumerate(args.vars)})

        

    tf.reset_default_graph()

outputs_ensemble = {base_name[I]:np.mean([all_outputs[J][base_name[I]] for J,j in enumerate(models)],0) for I,i in enumerate(args.vars)}
thresholds = {'MCC':0.426,'Sw':0.084}

if args.quiet == 0:
    print('Writing files...')
for J,j in enumerate(tqdm.tqdm(base_name,file=tqdmfile)):
    with open(j+'.spotds','w') as f:
        f.write('# SPOT-Disorder2-Seq prediction for %s\n'%(args.vars[J]))
        f.write('# Threshold for maximizing MCC: %1.3f'%(thresholds['MCC']))
        f.write('# Threshold for maximizing Sw: %1.3f'%(thresholds['Sw']))
        f.write('Pos\tAA\tP(D)\tLab\n')
        if args.vars[J][-3:] == '.seq':
            seq = read_seq_file(args.vars[J])
        else:
            seq = read_fasta_file(args.vars[J])
        for K,k in enumerate(seq):
            label = 'D' if outputs_ensemble[j][K]>=thresholds['MCC'] else 'O'
            f.write('%i\t%s\t%1.4f\t%s\n'%(K+1,k,outputs_ensemble[j][K],label))
        
        


