def translate_aa(three_letter):
 trans = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
 return (trans[three_letter])


def job_submit_pssm(file_name, dir):
 if os.path.exists(dir+'/'+file_name+'.pssm'):
  return
 print("Build "+file_name+" PSSM...")
 global BLASTDATABASE
 global BLAST_NUM_THREADS
 print(BLASTDATABASE, BLAST_NUM_THREADS)
 cmd = "psiblast -query "+dir+'/'+file_name+'.fasta'+" -num_threads "+str(BLAST_NUM_THREADS)+" -db "+BLASTDATABASE+" -num_iterations 3 -out "+dir+'/'+file_name+".out -out_ascii_pssm "+dir+'/'+file_name+".pssm 2>/dev/null"
 out = getstatusoutput(cmd)
 if out[0] != 0:
     print('Error: pssm profile not generated')

def job_submit_spd33(file_name, dir):
 if os.path.exists(dir+'/'+file_name+'.spd33'):
  return
 print("Build "+file_name+" spd33...")
 cmd = './run_list_spd33.sh '+file_name+' '+dir
 out = getstatusoutput(cmd)
 if out[0] != 0:
     print('error spd33 file not generated')

def job_submit_spotdis(file_name, dir):
 if os.path.exists(dir+'/'+file_name+'.spotds'):
  return
 print("Build "+file_name+" Spot Disorder Single...")
 cmd = './run_list_spotd.sh '+dir+'/'+file_name
 out=getstatusoutput(cmd)
 if out[0]!=0:
     print('error spotd file not generated')
 #cmd = './run_spotdis_single.py --gpu 0 --batch_size 50 --quiet 1 --output_dir seq_train seq_train/%s' % file_name+'.fasta'


def job_submit_alphafold2(file_name, dir):
 if os.path.exists(dir+'/'+file_name+'.af2.pdb'):
  return
 print("Build "+file_name+" Alpha Fold2 predicted structure...")
 cmd = './run_list_colab_alphafold2.sh '+dir+'/'+file_name+'.fasta '+dir+'/ '+file_name
 out=getstatusoutput(cmd)
 if out[0]!=0:
     print('error AlphaFold2, file not generated')


def pssm_check(file_name, dir):
 file_name += '.pssm'
 #print(file_name)
 if os.path.exists(dir+'/'+file_name):
  print("Build "+file_name+" PSSM done...")
 else:
  print("Can't get the PSSM, please check the sequence")
  sys.exit()

def spd33_check(file_name, dir):
 file_name += '.spd33'
 #print(file_name)
 if os.path.exists(dir+'/'+file_name):
  print("Build "+file_name+" spd33 done...")
 else:
  print("Can't get the spd33, please check the sequence")
  sys.exit()

def spotdis_check(file_name, dir):
 file_name += '.spotds'
 #print(file_name)
 if os.path.exists(dir+'/'+file_name):
  print("Build "+file_name+" Spot Disorder Single done...")
 else:
  print("Can't get the Spot Disorder Single, please check the sequence")
  sys.exit()

def alphafold2_check(file_name, dir):
 file_name += '.af2.pdb'
 #print(file_name)
 #print(dir+'/PDB/'+file_name)
 if os.path.exists(dir+'/'+file_name):
  print("Build "+file_name+" AlhpaFold2 structure prediction done...")
 else:
  print("Can't get the predicted structure from AlphaFold2, please check the sequence")
  sys.exit()


def pharmcophore_count_diff(wild_pharm_count, mutant_pharm_count):
 """Calculate the difference between to pharmacophore counts"""
 atom_class = ('A', 'D', 'H')
 pharm_w = [0] * 3
 i = 0
 for item in atom_class:
  if item in wild_pharm_count:
   pharm_w[i] = wild_pharm_count.get(item)
  i += 1
 pharm_m = [0] * 3
 i = 0
 for item in atom_class:
  if item in mutant_pharm_count:
   pharm_m[i] = mutant_pharm_count.get(item)
  i += 1
 return [a_i - b_i for a_i, b_i in zip(pharm_w, pharm_m)]


def features(protein_chain,mutation_pos,wild_aa,mutation_aa,mutation_list, a, b, c, Temp, pH, dir):
 '''Extracting features'''
 global BLASTDATABASE
 global BLAST_NUM_THREADS
 global uniclust30
 BLASTDATABASE=a
 BLAST_NUM_THREADS=b
 uniclust30=c
 if len(wild_aa) == 3:
  wild_aa=translate_aa(wild_aa)
 if len(mutation_aa) ==3:
  mutation_aa=translate_aa(mutation_aa)
 label = []
 if Temp:
     label.append(Temp)
 else:
     label.append('25')
 if pH:
     label.append(pH)
 else:
     label.append('7')
 label.append(muttype_ala.muttype_ala(mutation_aa))
 label.append(net_volume.net_volume(wild_aa,mutation_aa))
 label.append(net_hydrophobicity.net_hydrophobicity(wild_aa,mutation_aa))
 label.append(flexibility.flexibility(wild_aa,mutation_aa))
 label.append(mutation_hydrophobicity.mutation_hydrophobicity(wild_aa,mutation_aa))
 label.append(mutation_polarity.mutation_polarity(wild_aa,mutation_aa))
 label.append(mutation_type.mutation_type(wild_aa,mutation_aa))
 label.append(mutation_size.mutation_size(wild_aa,mutation_aa))
 label.append(mutation_hbonds.mutation_hbonds(wild_aa,mutation_aa))
 label.append(mutation_chemical.mutation_chemical(wild_aa,mutation_aa))
 label = label+sequences_features.sequences_features(protein_chain, mutation_pos, wild_aa, dir)
 if mutation_list==0:
  job_submit_pssm(protein_chain, dir)
  job_submit_spd33(protein_chain, dir)
  spd33_check(protein_chain, dir)
  pssm_check(protein_chain, dir)
  job_submit_spotdis(protein_chain, dir)
  spotdis_check(protein_chain, dir)
  job_submit_alphafold2(protein_chain, dir)
  alphafold2_check(protein_chain, dir)
 label=label+mpssmscores.mpssmscores(protein_chain,mutation_pos, dir)
 label=label+psepssm.psepssm(protein_chain, dir)
 #BoostddG features
 label.append(get_dF.get_dF(protein_chain, mutation_pos, wild_aa, mutation_aa, dir))
 label = label+get_daaph7.get_daaph7(wild_aa, mutation_aa)
 label.append(cal_aa_score.cal_aa_score(wild_aa, mutation_aa))
 label = label+parse_spd33.parse_spd33(protein_chain, mutation_pos, wild_aa, dir)
 #label = label+parse_spd33(protein_chain, mutation_pos, wild_aa, dir)
 label.append(get_spotd.get_spotd(protein_chain, mutation_pos+1, dir))
 #iFeature features
 label = label+z_scale.z_scale(wild_aa, mutation_aa)
 label.append(modif_gaac.modif_gaac(wild_aa, mutation_aa))
 #DDGun Features
 label= label+ddgun_features(dir+'/'+protein_chain+'.fasta', wild_aa, mutation_pos, mutation_aa, c, dir)
 label= label+mCSM_features.dssp_features(dir+'/'+protein_chain+'.af2.pdb', 'A', int(mutation_pos))
 pharm_sign, pharm_count = mCSM_features.pahrmaco_sign(dir + '/', protein_chain + '.af2', 'A', wild_aa,
                                                       mutation_pos)
 label = label+pharm_sign
 #print(pharm_sign, pharm_count)

 return label, pharm_count

import sys, getopt
import os
import numpy as np
import xgboost as xgb
from itertools import islice
import re
import math
import joblib
import subprocess
import shlex
from features import *
from ddgun_seq import ddgun_features
from subprocess import getstatusoutput
