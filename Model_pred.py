import ml


def read_mutfile(file_name):
    '''Reading mutation file'''
    mut_list = []
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            mut_list.append(line.strip().split(' '))
    file.close()
    return mut_list

def seqfile_prep(uniq_id, seq_header, sequence, Gfile_dir):
    '''Preparing fasta file for feature extraction.'''
    f = open(Gfile_dir + '/' + uniq_id + '.fasta', "w+")
    f.write(">%s.fasta \n%s" % (uniq_id, sequence))
    f.close()
    f = open(Gfile_dir + '/' + uniq_id, "w+")
    f.write("%s\n%s" % (seq_header, sequence))
    f.close()

def read_seq(fasta_seq, Gfile_dir):
    '''Read sequence file, make it ready for feature extraction, and generate a unique ID.
    Input: Fasta file
    Outpt: A unique ID'''
    seq_header = ''
    sequence = ''
    seq = open(fasta_seq, 'r')
    for line in seq:
        if line[0] == ">":
            seq_header = line.strip()
        else:
            sequence = sequence + line.strip()
    sequence = str(sequence)
    seq_header = str(seq_header)
    from hashlib import blake2b
    h = blake2b(digest_size=4)
    h.update(str(sequence).encode('utf-8'))
    uniq_id = h.hexdigest()
    uniq_id = str(uniq_id)
    seqfile_prep(uniq_id, seq_header ,sequence, Gfile_dir)
    return uniq_id


#convert one-letter amino acid to three-letter format
def translate_1aa3(one_letter):
    trans = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS',
             'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
             'Y': 'TYR', 'V': 'VAL'}
    return (trans[one_letter])


def get_options():
    '''Argument parsing'''
    parser = argparse.ArgumentParser(usage='Invalid arguments.', description='Program for predicting the change in stability (∆∆G(kcal/mol)) upon single point missense mutation.')
    parser.add_argument('-file','--file', required=True, help='Input sequence file in FASTA format')
    parser.add_argument('-mutation','--mutation', nargs='+', type=str, help='Inline mutation (example: A 32 G)')
    parser.add_argument('-mutlist','--mutlist','--ml', '--mutation-list', action='store', type=str, dest='ml', help='Mutation-list (as provided in the example directory).')
    parser.add_argument('-outdir',"--outdir", "--out-dir", action="store",type=str, dest="outdir", help="Output directory")
    parser.add_argument("-out-file", "--out-file", action="store", type=str, dest="outfile", help="Output result file")
    args = parser.parse_args()
    if not os.path.isfile(args.file):
        print('Error: Incorrect sequence file')
        sys.exit(1)
    seq_id = read_seq(args.file, 'Gen_Files')
    print('A unique id for the sequence:', seq_id)
    mut_file =None
    if args.ml:
        mut_data = read_mutfile(args.ml)
    else:
        mut_data =[]
        mut_data.append(args.mutation)
    outdir = os.getcwd()
    outfile = 'Result'
    if args.outdir: outdir = args.outdir
    if args.outfile: outfile = args.outfile
    return seq_id, mut_data, outdir, outfile


def feature_label_gen():
    '''The following feature list is used.'''
    # features
    feature_label = []
    feature_label.append('Temp')
    feature_label.append('pH')
    feature_label.append('Ala/NonAla')
    feature_label.append('NET_VOL')
    feature_label.append('NET_HPHO')
    feature_label.append('FLEX')
    feature_label.append('MUT_HPHO')
    feature_label.append('MUT_POL')
    feature_label.append('MUT_TYPE')
    feature_label.append('MUT_SIZE')
    feature_label.append('MUT_HHbonds')
    feature_label.append('MUT_CHEM')
    for i in range(0, 11):
        feature_label.append('RES_Neighb' + str(i))
    for i in range(0, 140):
        feature_label.append('PSSM' + str(i))
    for i in range(0, 160):
        feature_label.append('psePSSM' + str(i))
    feature_label.append('dF')
    for i in range(0, 7):
        feature_label.append('daaph' + str(i))
    feature_label.append('aa_score')
    for i in range(0, 4):
        feature_label.append('spd33' + str(i))
    feature_label.append('spotd')
    for i in range(0, 5):
        feature_label.append('z_scale' + str(i))
    feature_label.append('Modif_GAAC')
    for i in range(0, 3):
        feature_label.append('DDGun_feature' + str(i))
    feature_label.append('RSA')
    feature_label.append('phi')
    feature_label.append('psi')
    for i in range(0, 12):
        feature_label.append('mCSM_sign(' + str(i) + ')')
    for i in range(0, 3):
        feature_label.append('pharmaco_diff' + str(i))
    return feature_label

def main(a, b, c):

 #parse arguments
 seq_id, mut_data, outdir, outfile = get_options()
 #feature labels
 feature_label = feature_label_gen()

 params = {'colsample_bytree': 0.5, 'learning_rate': 0.05, 'max_depth': 9, 'alpha': 1.2, 'lambda': 2.0, 'gamma': 0.1,
           'min_child_weight': 6, 'min_samples_leaf': 1, 'min_samples_split': 2, 'subsample': 1,
           'max_features': 0.18}
 #####optimizing max_depth and min_child_weight####
 '''gridsearch_params = [
     (max_depth, min_child_weight)
     for max_depth in range(7, 11)
     for min_child_weight in range(5, 9)
 ]
 min_mae = float("Inf")
 best_params = None
 for max_depth, min_child_weight in gridsearch_params:
     print ("CV with max_depth={}, min_child_weight={}".format(
                             max_depth,
                             min_child_weight))
     # Update our parameters
     params['max_depth'] = max_depth
     params['min_child_weight'] = min_child_weight
     cv_results = xgb.cv(
         params,
         dtrain=data_dmatrix,
         num_boost_round=num_boost_round,
         seed=42,
         nfold=5,
         metrics={'rmse'},
         early_stopping_rounds=50
     )
     mean_mae = cv_results['test-rmse-mean'].min()
     boost_rounds = cv_results['test-rmse-mean'].argmin()
     print("\tRMSE {} for {} rounds".format(mean_mae, boost_rounds))
     if mean_mae < min_mae:
         min_mae = mean_mae
         best_params = (max_depth, min_child_weight)
 print("Best params: {}, {}, RMSE: {}".format(best_params[0], best_params[1], min_mae))'''
 ####optimizing colsample_bytree and subsample#####
 '''gridsearch_params = [
    (subsample, colsample)
    for subsample in [i/10. for i in range(8,11)]
    for colsample in [i/10. for i in range(3,8)]
 ]

 # Define initial best params and RMSE
 min_rmse = float("Inf")
 best_params = None
 for subsample, colsample in reversed(gridsearch_params):
     print("CV with subsample={}, colsample_bytree={}".format(
                             subsample,
                             colsample))

     #Update our parameters
     params['subsample'] = subsample
     params['colsample_bytree'] = colsample
     # Run CV
     cv_results = xgb.cv(
         params,
         data_dmatrix,
         num_boost_round=num_boost_round,
         seed=42,
         nfold=5,
         metrics={'rmse'},
         early_stopping_rounds=800
     )
     # Update best MAE
     mean_rmse = cv_results['test-rmse-mean'].min()
     boost_rounds = cv_results['test-rmse-mean'].argmin()
     print("\tRMSE {} for {} rounds".format(mean_rmse, boost_rounds))
     if mean_rmse < min_rmse:
         min_rmse = mean_rmse
         best_params = (subsample, colsample)
 print("Best params: {}, RMSE: {}".format(best_params, min_rmse))'''

 # Learning rate optimization
 '''min_rmse = float("Inf")
 best_params = None
 for lrt in [0.02, 0.03, 0.04, .07]:
     print("CV with learning rate={}".format(lrt))
     # We update our parameters
     params['learning_rate'] = lrt
     # Run and time CV
     cv_results = xgb.cv(
         params,
         dtrain=data_dmatrix,
         num_boost_round=num_boost_round,
         seed=42,
         nfold=5,
         metrics=['rmse'],
         early_stopping_rounds=100
     )
     # Update best score
     mean_rmse = cv_results['test-rmse-mean'].min()
     boost_rounds = cv_results['test-rmse-mean'].argmin()
     print("\tRMSE {} for {} rounds\n".format(mean_rmse, boost_rounds))
     if mean_rmse < min_rmse:
         min_rmse = mean_rmse
         best_params = lrt
 print("Best params: {}, RMSE: {}".format(best_params, min_rmse))'''
 ###Alhpa and lambda optimization
 '''gridsearch_params = [
      (reg_alpha, reg_lambda)
      for reg_alpha in [i/10. for i in range(10, 15)]
      for reg_lambda in [i/10. for i in range(15, 21)]
  ]
 min_rmse = float("Inf")
 best_params = None
 for reg_alpha, reg_lambda in gridsearch_params:
     print ("CV with reg_alpha={}, reg_lambda={}".format(reg_alpha,reg_lambda))
     #Update our parameters
     params['alpha'] = reg_alpha
     params['lambda'] = reg_lambda
     cv_results = xgb.cv(
         params,
         dtrain=data_dmatrix,
         num_boost_round=num_boost_round,
         seed=42,
         nfold=5,
         metrics={'rmse'},
         early_stopping_rounds=750)
     mean_rmse = cv_results['test-rmse-mean'].min()
     boost_rounds = cv_results['test-rmse-mean'].argmin()
     print("\tRMSE {} for {} rounds".format(mean_rmse, boost_rounds))
     if mean_rmse < min_rmse:
         min_rmse = mean_rmse
         best_params = (reg_alpha, reg_lambda)
 print("Best params: {}, {}, RMSE: {}".format(best_params[0], best_params[1], min_rmse))'''
 #optimizing min_sample_split and min_sample_leaf
 '''gridsearch_params = [
     (min_leaf, min_split)
     for min_leaf in range(1, 10)
     for min_split in range(2, 10)
 ]
 min_rmse = float("Inf")
 best_params = None
 for min_leaf, min_split in gridsearch_params:
     print ("CV with min_leaf={}, min_split={}".format(
         min_leaf,
         min_split))
     # Update our parameters
     params['min_samples_leaf'] = min_leaf
     params['min_samples_split'] = min_split
     cv_results = xgb.cv(
         params,
         dtrain=data_dmatrix,
         num_boost_round=num_boost_round,
         seed=42,
         nfold=5,
         metrics={'rmse'},
         early_stopping_rounds=750
     )
     mean_rmse = cv_results['test-rmse-mean'].min()
     boost_rounds = cv_results['test-rmse-mean'].argmin()
     print("\tRMSE {} for {} rounds".format(mean_rmse, boost_rounds))
     if mean_rmse < min_rmse:
         min_rmse = mean_rmse
         best_params = (min_leaf, min_split)
 print("Best params: {}, {}, RMSE: {}".format(best_params[0], best_params[1], min_rmse))'''
 #Gamma optimization
 '''min_rmse = float("Inf")
 best_params = None
 for gm in [0.1, 0.2, 0.3, 0.4]:
     print("CV with gamma ={}".format(gm))
     #We update our parameters
     params['gamma'] = gm
     #Run and time CV
     cv_results = xgb.cv(
          params,
          dtrain=data_dmatrix,
          num_boost_round=num_boost_round,
          seed=42,
          nfold=5,
          metrics=['rmse'],
          early_stopping_rounds=750
     )
     #Update best score
     mean_rmse = cv_results['test-rmse-mean'].min()
     boost_rounds = cv_results['test-rmse-mean'].argmin()
     print("\tRMSE {} for {} rounds\n".format(mean_rmse, boost_rounds))
     if mean_rmse < min_rmse:
         min_rmse = mean_rmse
         best_params = gm
 print("Best params: {}, RMSE: {}".format(best_params, min_rmse))'''
 #Final Cross Validation
 '''cv_results = xgb.cv(dtrain=data_dmatrix, params=params, nfold=5,
                     num_boost_round=300, as_pandas=True, seed=42, early_stopping_rounds=150, metrics=['rmse'])
 print(cv_results.head)
 print((cv_results["test-rmse-mean"]).idxmin(), (cv_results["test-rmse-mean"]).min())'''
 #Visualizing Tree and feature Importance, SAVE MODEL
 #xg_reg = xgb.train(params=params, dtrain=data_dmatrix, num_boost_round=300)
 #xgb.plot_tree(xg_reg, num_trees=1)
 ##plt.show()
 #xgb.plot_importance(xg_reg)
 #plt.rcParams['figure.figsize'] = [10, 10]
 #plt.show()
 #joblib.dump(xg_reg, 'S2648reverse')
 #xg_reg.save_model("S2648reverse_spd_spotd.model")
 '''testdata_list = [['CAGI5_PTENboost','CAGI5_TPMT_PTEN', 'CAGI5_TPMT_PTEN'],
                  ['CAGI5_TPMTboost', 'CAGI5_TPMT_PTEN', 'CAGI5_TPMT_PTEN'],
                  ['CAGI5_TPMT_PTENboost', 'CAGI5_TPMT_PTEN', 'CAGI5_TPMT_PTEN'],
                  ['Ssym_direct', 'Ssympdbsorted_direct', 'Ssym_direct'],
                  ['Ssym_inverse', 'Ssympdbsorted_inverse', 'Ssym_inverse']]
 for tst in testdata_list:'''
 gen_file_dir = 'Gen_Files' #storage of internal files
 test_model(a, b, c, seq_id, mut_data, gen_file_dir, feature_label, outdir, outfile)


def test_model(a, b, c, uniq_id, mut_data, dir, feature_label, outdir, outfile):
 #Testing Model
 pdbidc_last = 'NA'
 features_list = []
 mut_list = []
 for r in mut_data:
    mut_list.append(r[0]+' '+r[1]+' '+r[2])
    mut_seqpos = int(r[1])
    wild_res = r[0]
    mut_res = r[2]
    mut_info = translate_1aa3(wild_res) + '-' + str(mut_seqpos) + '-' + translate_1aa3(mut_res)
    if uniq_id != pdbidc_last:
        pdbidc_last = uniq_id
        wildtype_features, wild_pharm_count = Feature_extract.features(uniq_id, mut_seqpos, wild_res, mut_res, 0, a, b, c, '25', '7', dir)
        #features_list.append(Feature_extract.features(pdbidc, mut_seqpos, mut_res, wild_res, 0, a, b, c, r[8], r[9], neg_ddg))
    else:
        wildtype_features, wild_pharm_count = Feature_extract.features(uniq_id, mut_seqpos, wild_res, mut_res, 1, a, b, c, '25', '7', dir)
        #features_list.append(Feature_extract.features(pdbidc, mut_seqpos, mut_res, wild_res, 1, a, b, c, r[8], r[9], neg_ddg))
    pdbidc_PosMutRes = uniq_id + '_' + str(mut_seqpos) + mut_res
    if os.path.exists(dir + '/' + pdbidc_PosMutRes + '.af2.pdb'):
        print('Mutated pdb already existed')
    else:
        pdb_mutate(dir + '/', uniq_id, 'A', mut_info, mut_seqpos, mut_res)
    #muttype_features, mutant_pharm_count = Feature_extract.features(uniq_id, mut_seqpos, mut_res, wild_res, 1, a, b, c, '25', '7', dir)
    muttype_pharm_sign, mutant_pharm_count = mCSM_features.pahrmaco_sign(dir + '/', uniq_id + '.af2', 'A', wild_res, mut_seqpos)
    #print(mutant_pharm_count)
    wildtype_features = wildtype_features + Feature_extract.pharmcophore_count_diff(wild_pharm_count, mutant_pharm_count)
    features_list.append(wildtype_features)

 #preparing data for testing
 features_list=list(np.array(features_list, dtype=np.float32))
 data = pd.DataFrame(features_list)
 data.columns = feature_label
 test_dmatrix = xgb.DMatrix(data)
 #prediction with weighted average ensemble model
 Xfeatures = np.array(data)
 ypred = ml.ens_predict(Xfeatures)
 #ypred = xg_reg.predict(test_dmatrix)
 '''from sklearn.metrics import mean_squared_error
 print('PCC', np.corrcoef(ytest, ypred)[0, 1])
 print('RMSE', np.sqrt(((ypred- ytest) ** 2).mean()))
 print('MSE', mean_squared_error(ytest, ypred))
 print('MAE', np.mean(np.absolute(ypred-ytest)))'''
 pred_ddG = list(ypred)
 mut_pred_ddG = []
 for i in range(len(mut_list)):
     mut_pred_ddG.append([mut_list[i],pred_ddG[i]])
 '''f = open(dir+'/'+uniq_id+'_result.txt', "w+")
 print('#Mut\t∆∆G(kcal/mol)')
 f.write('Mut\t∆∆G(kcal/mol)\n')
 for i in range(len(mut_list)):
     f.write("%s\t%s\n" %(mut_list[i], pred_ddG[i]))
 f.close()'''
 f = open(outdir + '/' + outfile + '.txt', "w+")
 print('#Mut\t∆∆G(kcal/mol)')
 f.write('Mut\t∆∆G(kcal/mol)\n')
 for i in range(len(mut_list)):
     f.write("%s\t%s\n" % (mut_list[i], pred_ddG[i]))
     print(mut_list[i], '\t', pred_ddG[i])
 f.close()
 #Feature Importance using xgboost
 #feature_imp = xg_reg.get_score(importance_type='gain')
 #print (feature_imp)
 '''for key, value in feature_imp.items():
     print (key,"\t",value)'''
 #Feature Importance using Mutual Information metric
 '''from sklearn.feature_selection import SelectKBest
 import matplotlib.pyplot as plt
 from sklearn.feature_selection import mutual_info_regression
 # training and feature selection
 MIf_selector = SelectKBest(score_func=mutual_info_regression, k='all').fit(f_datatrain, f_datalabel)
 #fitting a model
 #MIf_selector.fit =(f_datatrain, f_datalabel)
 #transform train and test input
 X_MIfs = MIf_selector.transform(f_datatrain)
 #Xtest_MIfs = MIf_selector.transform(Xtest)
 # Plot the scores for the features
 plt.bar([i for i in range(len(MIf_selector.scores_))], MIf_selector.scores_)
 plt.xlabel("feature index")
 plt.ylabel("Estimated MI value")
 plt.show()
 print('####Feature Importance using Mutual Information metric####')
 for i in MIf_selector.scores_:
     print(i)
#https://towardsdatascience.com/how-to-perform-feature-selection-for-regression-problems-c928e527bbfa
 ##Feature selection using Correlation metric##
 from sklearn.feature_selection import f_regression
 f_selector = SelectKBest(score_func=f_regression, k='all').fit(f_datatrain, f_datalabel)
 X_MIfs = f_selector.transform(f_datatrain)
 print('####Feature Importance using Correlation Metric####')
 for i in f_selector.scores_:
     print(i)
 ##Feature selection uisng MRMR##
 from sklearn.feature_selection import f_regression
 # inputs:
 #    X: pandas.DataFrame, features
 #    y: pandas.Series, target variable
 #    K: number of features to select
 X = pd.DataFrame(X)
 y = pd.Series(y)
 # compute F-statistics and initialize correlation matrix
 F = pd.Series(f_regression(X, y)[0], index = X.columns)
 corr = pd.DataFrame(.00001, index = X.columns, columns = X.columns)

 # initialize list of selected features and list of excluded features
 selected = []
 not_selected = X.columns.to_list()
 K = len(X.columns)
 K= 330
 # repeat K times
 for i in range(K):

     # compute (absolute) correlations between the last selected feature and all the (currently) excluded features
     if i > 0:
         last_selected = selected[-1]
         corr.loc[not_selected, last_selected] = X[not_selected].corrwith(X[last_selected]).abs().clip(.00001)

    # compute FCQ score for all the (currently) excluded features (this is Formula 2)
     score = F.loc[not_selected] / corr.loc[not_selected, selected].mean(axis = 1).fillna(.00001)

    # find best feature, add it to selected and remove it from not_selected
     best = score.index[score.argmax()]
     selected.append(best)
     not_selected.remove(best)

 print('Selected Features:', selected)
 print('Features ignored:', not_selected)'''

import argparse
import sys, getopt
import os
import numpy as np
import xgboost as xgb
import json
import pandas as pd
import joblib
import Feature_extract
from mutate_pdb import pdb_mutate
from features import mCSM_features
import pickle
import matplotlib.pyplot as plt

