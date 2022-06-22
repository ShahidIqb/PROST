#!/bin/bash

###...activate colab AlphaFold2 environment for colabfold_batch running...###
source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate /Users/siqb0003/PycharmProjects/colab_alphafold_local/colabfold_batch/colabfold-conda

mkdir $2temp_alphafold2
colabfold_batch --templates --num-recycle 3 --use-gpu-relax --cpu --num-models 1 $1 $2temp_alphafold2
mv $2temp_alphafold2/*.pdb $2$3.af2.pdb
rm -rf $2temp_alphafold2
