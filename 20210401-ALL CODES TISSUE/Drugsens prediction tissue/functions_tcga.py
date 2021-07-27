
import os
import numpy as np
import seaborn as sns
from sklearn import metrics

from my_load_data import * 
from my_util import *
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math



def my_tsne_RNA(X,y_map,the_label,the_color,title_name = None,flag_fpkm_tpm=False,bool_pca=True,pca_dim=30,logg=True,bool_save=False,save_path=''):
    
    if flag_fpkm_tpm:
        X1 = fpkm_to_tpm(X)
    else:
        X1=X
    if logg:
        X1 = np.log(X1+1) 

    X_normalized = X1

    # pca
    print(X_normalized.shape)
    if bool_pca and X1.shape[1]>pca_dim:
        pca=PCA(n_components=min(pca_dim, X1.shape[0]))
        X_pca = pca.fit_transform(X_normalized)
        print('sum(pca.explained_variance_ratio_)',sum(pca.explained_variance_ratio_))
    else:
        X_pca = X_normalized
    print(X_pca.shape)

    # tsne
    X_tsne = TSNE(n_components=2).fit_transform(X_pca)
    print(X_tsne.shape)

    
    for i in range(len(the_label)):
        plt.scatter(X_tsne[y_map==i,0], X_tsne[y_map==i,1],c=the_color[i], label = the_label[i])
    sns.despine()
    plt.legend(fontsize='small')
    plt.xlabel('t-SNE dim1')
    plt.ylabel('t-SNE dim2')
    if title_name!=None:
        plt.title(title_name)
    if bool_save and len(save_path)>0:
        plt.savefig(save_path+'.eps',format='eps')
        plt.savefig(save_path+'.png',dpi=600)
    plt.show()
    
    return 0

def my_pca_RNA(X,y_map,the_label,the_color,title_name = None,flag_fpkm_tpm=False,pca_dim=2,logg=True,bool_save=False,save_path=''):
    if flag_fpkm_tpm:
        X1 = fpkm_to_tpm(X)
    else:
        X1=X
    if logg:
        X1 = np.log(X1+1) 

    X_normalized = X1

    # pca
    print(X_normalized.shape)
    pca=PCA(n_components=min(pca_dim, X1.shape[0]))
    X_pca = pca.fit_transform(X_normalized)
    print('sum(pca.explained_variance_ratio_)',sum(pca.explained_variance_ratio_))
    print(X_pca.shape)
    
    for i in range(len(the_label)):
        plt.scatter(X_pca[y_map==i,0], X_pca[y_map==i,1],c=the_color[i], label = the_label[i])
    sns.despine()
    plt.legend(fontsize='small')
    plt.xlabel('PCA dim1')
    plt.ylabel('PCA dim2')
    if title_name!=None:
        plt.title(title_name)
    if bool_save and len(save_path)>0:
        plt.savefig(save_path+'.eps',format='eps',dpi=600)
        plt.savefig(save_path+'.png',dpi=600)
    plt.show()
    
    return 0



# knn_cv_predict
def get_knn_cv_predict_results(X, Y, C_neighbor_k, C_cv_fold, flag_display=False, drug_name='', flag_tpm_normalize=False):
    
    # basic checks
    if Y.shape[0] < C_cv_fold or X.shape[1]<1:  # not enough data
        return math.nan, math.nan, math.nan
    if len(np.unique(Y)) < 2:   # Y is not binarized
        return math.nan, math.nan, math.nan
    if C_neighbor_k > Y.shape[0]//2:   # C_neighbor_k
        C_neighbor_k = Y.shape[0]//2
    
    # TPM normalization
    if flag_tpm_normalize:
        X_sum = np.sum(X, axis=1) / 1e6
        X_sum[X_sum<1e-6] = 1
        for i in range(X.shape[0]):
            X[i, :] = X[i, :] / X_sum[i]
        X=np.log(X)
    
    # knn cross validation
    model_knn = KNeighborsClassifier(n_neighbors = C_neighbor_k, weights='distance', algorithm = 'auto')
    y_pred = cross_val_predict(model_knn, X, Y, cv = C_cv_fold)
    y_pred_proba = cross_val_predict(model_knn, X, Y, cv = C_cv_fold,method = 'predict_proba')[:,1]
    
    mse_score = np.mean((Y - y_pred)**2)
    
    acc_score = metrics.accuracy_score(Y,y_pred)
    
    auc_score = math.nan
    if len(np.unique(Y)) == 2:
        auc_score = metrics.roc_auc_score(Y, y_pred_proba)
            
    if flag_display:
        print('len_X',len(X),'C_neighbor_k',C_neighbor_k,'Accuracy=', acc_score, ';AUC=', auc_score)
    
    return y_pred, auc_score, mse_score

###################################################################################

# KNN enumeration
def get_knn_para_by_bruteforce_TCGA(X_idx_raw, X_raw, Y_raw, Y_drug_idx, C_cv_fold, flag_tpm_normalize=False,C_selected_para_max=30, C_selected_para_min=10, flag_second_best=False, second_best_th=0.05, bool_normalize = True):
    best_best_mean_auc_score = 0
    best_C_selected_para_cnt = 1
    best_best_C_neighbor_k = 1
    
    best_mean_auc_score_list = []
    best_C_neighbor_k_list = []
    
    
    if C_selected_para_max < C_selected_para_min:
        C_selected_para_max = C_selected_para_min
    C_selected_para_cnt_list = [i+C_selected_para_min for i in                                 range(min(C_selected_para_max+1-C_selected_para_min,                                           len(X_idx_raw)+1-C_selected_para_min))]
    for C_selected_para_cnt in C_selected_para_cnt_list:  # gene numbers

        X_idx = X_idx_raw[:C_selected_para_cnt]
        
        C_neighbor_k_max = 30
        C_neighbor_k_list = [i+1 for i in range(C_neighbor_k_max)]
        
        best_mean_auc_score = 0
        best_C_neighbor_k = 1
        
        for C_neighbor_k in C_neighbor_k_list:   # parameter k
        
            auc_score_list = []
            mse_score_list = []

            
            for i in range(len(Y_drug_idx)):  # loop of every drugs
                        
                drug_idx = Y_drug_idx[i]

                # load X,Y
                X = X_raw[:, X_idx]
                Y = Y_raw[drug_idx, :]

                # normalize Y
                if bool_normalize:
                    Y = my_normalize(Y)

                ## cross validate predict
                y_pred, auc_score, mse_score = get_knn_cv_predict_results(X, Y, C_neighbor_k, C_cv_fold, False, '', flag_tpm_normalize)
                if not math.isnan(auc_score):
                    auc_score_list.append(auc_score)
                    mse_score_list.append(mse_score)
            
            cnt_expected_results = len(Y_drug_idx)
            
            # record the current best parameters
            if len(auc_score_list)>=cnt_expected_results*0.8 and best_mean_auc_score < min(auc_score_list):
                best_mean_auc_score = min(auc_score_list)
                best_C_neighbor_k = C_neighbor_k
            
        
        print('C_selected_para_cnt=%d, best_C_neighbor_k=%d, AUC=%f'%(C_selected_para_cnt,best_C_neighbor_k, best_mean_auc_score))

        # record the best parameters
        best_mean_auc_score_list.append(best_mean_auc_score)
        best_C_neighbor_k_list.append(best_C_neighbor_k)
        if best_best_mean_auc_score < best_mean_auc_score:
            best_best_mean_auc_score = best_mean_auc_score
            best_best_C_neighbor_k = best_C_neighbor_k
            best_C_selected_para_cnt = C_selected_para_cnt
    
    print('[best] C_selected_para_cnt=%d, best_C_neighbor_k=%d, AUC=%f'%(best_C_selected_para_cnt, best_best_C_neighbor_k, best_best_mean_auc_score))
    
    
    # Suboptimal decision
    if flag_second_best:
        for i in range(len(best_mean_auc_score_list)):
            if best_mean_auc_score_list[i] > best_best_mean_auc_score - second_best_th :
                best_C_selected_para_cnt = i + C_selected_para_min
                best_best_C_neighbor_k = best_C_neighbor_k_list[i]
                best_best_mean_auc_score = best_mean_auc_score_list[i]
                break
        print('[second best] C_selected_para_cnt=%d, best_C_neighbor_k=%d, AUC=%f'%(best_C_selected_para_cnt, best_best_C_neighbor_k, best_best_mean_auc_score))

    return best_C_selected_para_cnt, best_best_C_neighbor_k

