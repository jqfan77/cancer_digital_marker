# KNN + cross validation
# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn

import numpy as np
import math
import matplotlib.pyplot as plt
import random

from sklearn import metrics
from sklearn.model_selection import cross_val_predict
from sklearn.neighbors import KNeighborsRegressor

from my_util import *
from my_filter import *


###################################################################################

# KNN + cross validation
def get_knn_cv_predict_results(X, Y, C_neighbor_k, C_cv_fold, C_drug_sens_threshold, flag_display=False, drug_name='', flag_tpm_normalize=False):
    
    # check data
    if Y.shape[0] < C_cv_fold or X.shape[1]<1:  # not enough data points
        if flag_display:
            return math.nan, math.nan, math.nan, [], []
        else:
            return math.nan, math.nan, math.nan
    if len(np.unique(Y > C_drug_sens_threshold)) < 2:   # Y is not binarized
        if flag_display:
            return math.nan, math.nan, math.nan, [], []
        else:
            return math.nan, math.nan, math.nan
    if C_neighbor_k > Y.shape[0]//2:   # C_neighbor_k cannot be too large
        C_neighbor_k = Y.shape[0]//2
    
    # TPM normalize
    if flag_tpm_normalize:
        X_sum = np.sum(X, axis=1) / 1e6
        X_sum[X_sum<1e-6] = 1
        for i in range(X.shape[0]):
            X[i, :] = X[i, :] / X_sum[i]
    
    # create KNN model, and cross-validate
    model_knn = KNeighborsRegressor(n_neighbors = C_neighbor_k, weights = 'uniform', algorithm = 'auto')
    y_pred = cross_val_predict(model_knn, X, Y, cv = C_cv_fold)
    
    mse_score = np.mean((Y - y_pred)**2)
    
    auc_score = math.nan
    if len(np.unique(Y > C_drug_sens_threshold)) == 2:
        auc_score = metrics.roc_auc_score(Y > C_drug_sens_threshold, y_pred)
        if flag_display:
            print('AUC=', auc_score, '; RMSE=', mse_score**0.5)
    
    if flag_display:
        fpr, tpr, thresholds = metrics.roc_curve(Y < C_drug_sens_threshold, y_pred, pos_label=0)
        plt.plot(fpr,tpr)
        plt.title('%s gene_cnt=%d k=%d AUC=%.2f'%(drug_name, X.shape[1], C_neighbor_k, auc_score))
        plt.show()
    
    if flag_display:
        return y_pred, auc_score, mse_score, fpr, tpr
    else:
        return y_pred, auc_score, mse_score


###################################################################################

# Enumurate for good performance (AUC) on cross-validation 
def get_knn_para_by_bruteforce(X_idx_raw, X_raw, Y_raw, Y_drug_idx, \
                               C_cv_fold, C_drug_sens_threshold, \
                               dataset_test_ids, dataset_id, tissue_id_list, \
                               tissue_selected_id_list_list, \
                               flag_tpm_normalize=False, \
                               C_selected_para_max=30, C_selected_para_min=10,\
                               flag_second_best=False, second_best_th=0.05):
    best_best_mean_auc_score = 0
    best_C_selected_para_cnt = 1
    best_best_C_neighbor_k = 1
    
    best_mean_auc_score_list = []
    best_C_neighbor_k_list = []
    
    
    if C_selected_para_max < C_selected_para_min:
        C_selected_para_max = C_selected_para_min
    C_selected_para_cnt_list = [i+C_selected_para_min for i in \
                                range(min(C_selected_para_max+1-C_selected_para_min, \
                                          len(X_idx_raw)+1-C_selected_para_min))]
    for C_selected_para_cnt in C_selected_para_cnt_list:  # first `C_selected_para_cnt` will be used

        X_idx = X_idx_raw[:C_selected_para_cnt]
        
        C_neighbor_k_max = 30
        C_neighbor_k_list = [i+1 for i in range(C_neighbor_k_max)]
        
        best_mean_auc_score = 0
        best_C_neighbor_k = 1
        
        for C_neighbor_k in C_neighbor_k_list:   # `C_neighbor_k` for KNN
        
            auc_score_list = []
            mse_score_list = []

            for select_dataset in dataset_test_ids:  # specify dataset
                for tissue_selected_id_list in tissue_selected_id_list_list: # specify tissue
                    for i in range(len(Y_drug_idx)):  # specify drug
                        
                        drug_idx = Y_drug_idx[i]

                        # load X,Y
                        X = X_raw[:, X_idx]
                        X = np.log(1+X)
                        Y = Y_raw[drug_idx, :]
                        
                        # select cell lines with specified dataset and tissue
                        X_cell_idx = my_data_selection_2(dataset_id, select_dataset, tissue_id_list, tissue_selected_id_list, False)
                        X_cell_idx = X_cell_idx[~np.isnan(Y[X_cell_idx])]    # remove those without drug data
                        # get X,Y
                        X = X[X_cell_idx, :]
                        Y = Y[X_cell_idx]

                        # re-normalize Y to ensure mean=0 and std=1
                        Y = my_normalize(Y)

                        ## cross-validate 
                        y_pred, auc_score, mse_score = get_knn_cv_predict_results(X, Y, C_neighbor_k, C_cv_fold, C_drug_sens_threshold, False, '', flag_tpm_normalize)
                        if not math.isnan(auc_score):
                            auc_score_list.append(auc_score)
                            mse_score_list.append(mse_score)
                    #print('K=', C_neighbor_k, '; AUC=', np.mean(auc_score_list), '; RMSE=', np.mean(mse_score_list)**0.5)
            
            cnt_expected_results = len(dataset_test_ids) * len(tissue_selected_id_list_list) * len(Y_drug_idx)
            
            # record best K_KNN; compare min(auc_score_list);
            if len(auc_score_list)>=cnt_expected_results*0.8 and best_mean_auc_score < min(auc_score_list):
                best_mean_auc_score = min(auc_score_list)
                best_C_neighbor_k = C_neighbor_k
                
        print('.', end='')   # Running...
        
        # record best para_cnt and K_KNN
        best_mean_auc_score_list.append(best_mean_auc_score)
        best_C_neighbor_k_list.append(best_C_neighbor_k)
        if best_best_mean_auc_score < best_mean_auc_score:
            best_best_mean_auc_score = best_mean_auc_score
            best_best_C_neighbor_k = best_C_neighbor_k
            best_C_selected_para_cnt = C_selected_para_cnt
    
    print()
    print('[best1] C_selected_para_cnt=%d, best_C_neighbor_k=%d, AUC=%f'%(best_C_selected_para_cnt, best_best_C_neighbor_k, best_best_mean_auc_score))
    
    
    # find best2
    if flag_second_best:
        for i in range(len(best_mean_auc_score_list)):
            if best_mean_auc_score_list[i] > best_best_mean_auc_score - second_best_th :
                best_C_selected_para_cnt = i + C_selected_para_min
                best_best_C_neighbor_k = best_C_neighbor_k_list[i]
                best_best_mean_auc_score = best_mean_auc_score_list[i]
                break
        print('[best2] C_selected_para_cnt=%d, best_C_neighbor_k=%d, AUC=%f'%(best_C_selected_para_cnt, best_best_C_neighbor_k, best_best_mean_auc_score))

    return best_C_selected_para_cnt, best_best_C_neighbor_k

