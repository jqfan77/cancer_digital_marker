# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn

import math
import numpy as np
from scipy.stats import pearsonr


##############################################################################################
# select cell lines

# select cell lines with specified dataset, and tissue
# return index table `X_cell_idx`
def my_data_selection_2(dataset_id, select_dataset, tissue_id_list, tissue_selected_id_list, display_flag=True):
    X_tissue_id = tissue_id_list
    X_cell_idx = np.array([i for i in range(dataset_id.shape[0])])
    if select_dataset>0:
        X_tissue_id = X_tissue_id[dataset_id==select_dataset]
        X_cell_idx = X_cell_idx[dataset_id==select_dataset]
        if display_flag:
            print('selecting dataset',select_dataset,'  X_cell_idx:',np.shape(X_cell_idx))
    if len(tissue_selected_id_list)>0:
        tissue_sub_idx = np.zeros(X_tissue_id.shape[0], np.bool_)
        for tissue_id in tissue_selected_id_list:
            tissue_sub_idx = tissue_sub_idx | (X_tissue_id==tissue_id)
        X_cell_idx = X_cell_idx[tissue_sub_idx]
        if display_flag:
            print('selecting tissue; ', ' X_cell_idx:',np.shape(X_cell_idx))
    return X_cell_idx

# select cell lines with specified dataset, and tissue
# return gene expression `X`, drug zscore `Y`
def my_data_selection(X, Y, dataset_id, select_dataset, tissue_id_list, tissue_selected_id_list, display_flag=True):
    tissue_id_list_sub = tissue_id_list
    if select_dataset>0:
        X = X[dataset_id==select_dataset,:]
        Y = Y[dataset_id==select_dataset]
        tissue_id_list_sub=tissue_id_list_sub[dataset_id==select_dataset]
        if display_flag:
            print('selecting dataset',select_dataset,'  X:',np.shape(X), ', Y:',np.shape(Y))
    if len(tissue_selected_id_list)>0:
        tissue_sub_idx = np.zeros(Y.shape[0], np.bool_)
        for tissue_id in tissue_selected_id_list:
            tissue_sub_idx = tissue_sub_idx | (tissue_id_list_sub==tissue_id)
        X = X[tissue_sub_idx, :]
        Y = Y[tissue_sub_idx]
        if display_flag:
            print('selecting tissue; ',tissue_selected_id_list, ' X:',np.shape(X),', Y:',np.shape(Y))
    # remove those without drug data
    X = X[~np.isnan(Y),:]
    Y = Y[~np.isnan(Y)]
    if display_flag:
        print('remove cell without drug data: X:',np.shape(X),', Y:',np.shape(Y))
    return X,Y
    
    

##############################################################################################
# Gene filters (basic)

# return intersection of `idx_list_1` and `idx_list_2`
def get_repeating_idx_filter(idx_list_1, idx_list_2):
    Filter_idx = []
    '''for i in range(len(idx_list_1)):
        if idx_list_1[i] in idx_list_2:
            Filter_idx.append(i)'''
    for idx_2 in idx_list_2:
        if idx_2 in idx_list_1:
            Filter_idx.append(np.argwhere(idx_list_1==idx_2)[0][0])
    return np.array(Filter_idx)


# return `idx_list_1` - `idx_list_2`
def remove_repeating_idx_filter(idx_list_1, idx_list_2):
    Filter_idx = []
    for i in range(len(idx_list_1)):
        if not idx_list_1[i] in idx_list_2:
            Filter_idx.append(i)
    return np.array(Filter_idx)
    

##############################################################################################
# Gene filters (Linear)

# Calculate correlation efficient
# return genes with correlation efficient greater than `th`
def get_corr_filter_th(X, Y, th):
    Filter_idx = []
    for i in range(X.shape[1]):
        x_corr = pearsonr(X[:,i],Y)[0]
        if abs(x_corr)>th:
            Filter_idx.append(i)
    return Filter_idx

# Calculate correlation efficient
# return `n_para` genes
def get_corr_filter(X, Y, n_para):
    X_corr = np.zeros(X.shape[1])
    for i in range(X.shape[1]):
        X_corr[i] = pearsonr(X[:,i],Y)[0]
    X_corr_sorted_idx = sorted(range(X.shape[1]), key=lambda k: abs(X_corr[k]), reverse=True)
    X_corr_sorted = X_corr[X_corr_sorted_idx]
    
    return X_corr_sorted_idx[:n_para], X_corr_sorted[:n_para]

### calculate and sort standardized slope, negative correlation only
def get_slope_filter(X, Y):
    X = np.exp(X) - 1
    slope_list = np.zeros(X.shape[1])
    for i in range(X.shape[1]):
        X_sub = X[:, i]
        z = np.polyfit(Y, X_sub, 1)
        slope_list[i] = -z[0]
    slope_sorted_idx = sorted(range(len(slope_list)), key=lambda k:slope_list[k], reverse=True)
    slope_sorted = slope_list[slope_sorted_idx]
    # remove those lower than 0
    for i in range(len(slope_sorted)):
        if slope_sorted[i]<0:
            slope_sorted     = slope_sorted    [:i]
            slope_sorted_idx = slope_sorted_idx[:i]
            break
    return slope_sorted_idx, slope_sorted

# calculate and sort standardized slope
def get_std_slope_filter(X, Y):
    #X = np.exp(X) - 1
    slope_list = np.zeros(X.shape[1])
    for i in range(X.shape[1]):
        X_sub = X[:, i]
        # normalize X_sub
        X_sub = X_sub / np.mean(X_sub)
        
        z = np.polyfit(Y, X_sub, 1)
        slope_list[i] = z[0]
    slope_sorted_idx = sorted(range(len(slope_list)), key=lambda k:abs(slope_list[k]), reverse=True)
    slope_sorted = slope_list[slope_sorted_idx]
    # remove those lower than 0
    '''for i in range(len(slope_sorted)):
        if slope_sorted[i]<0:
            slope_sorted     = slope_sorted    [:i]
            slope_sorted_idx = slope_sorted_idx[:i]
            break'''
    return slope_sorted_idx, slope_sorted


##############################################################################################
# Gene filters (Dispersion)

# Fisher filter: compare the Fisher's discriminant of gene expression `X` between cell lines with the highest and lowest `top_bottom_ratio` zscore `Y`; it is one of the dispersion filter;
# return `n_para` genes with best scores
def get_topbottom_fisher_filter(X, Y, n_para, top_bottom_ratio):
    
    Y_sorted_idx = sorted(range(len(Y)), key=lambda k: Y[k])
    
    if top_bottom_ratio>0.5:
        top_bottom_ratio = 0.5
    Y_sorted_cnt = int(Y.size*top_bottom_ratio)

    Y_label = np.zeros(Y.shape)
    Y_label[Y_sorted_idx[:Y_sorted_cnt]] = 1       # bottom
    Y_label[Y_sorted_idx[-Y_sorted_cnt:]] = 2      # top  

    X_bottom = X[Y_label==1, :]
    X_top    = X[Y_label==2, :]
    X_bottom_mean = np.mean(X_bottom, axis=0)
    X_top_mean    = np.mean(X_top,    axis=0)
    X_bottom_std = np.std(X_bottom, axis=0, ddof=0)
    X_top_std    = np.std(X_top,    axis=0, ddof=0)
    X_top_bottom_score = (X_top_mean-X_bottom_mean) / (X_bottom_std+X_top_std+1e-5)

    X_top_bottom_score_sorted_idx = sorted(range(len(X_top_bottom_score)), key=lambda k: abs(X_top_bottom_score[k]), reverse=True)
    X_top_bottom_score_sorted = X_top_bottom_score[X_top_bottom_score_sorted_idx]
    
    return X_top_bottom_score_sorted_idx[:n_para], X_top_bottom_score_sorted[:n_para]

# Fisher filter: compare the Fisher's discriminant of gene expression `X` between cell lines with the highest `top_ratio` and lowest `bottom_ratio` zscore `Y`; it is one of the dispersion filter;
# return `n_para` genes with best scores
def get_top_bottom_fisher_filter(X, Y, n_para, bottom_ratio, top_ratio):
    Y_sorted_idx = sorted(range(len(Y)), key=lambda k: Y[k])
    
    Y_sorted_cnt_bottom = int(Y.size * bottom_ratio)
    Y_sorted_cnt_top    = int(Y.size * top_ratio)

    Y_label = np.zeros(Y.shape)
    Y_label[Y_sorted_idx[:Y_sorted_cnt_bottom]] = 1      # bottom
    Y_label[Y_sorted_idx[-Y_sorted_cnt_top:  ]] = 2      # top  

    X_bottom = X[Y_label==1, :]
    X_top    = X[Y_label==2, :]
    X_bottom_mean = np.mean(X_bottom, axis=0)
    X_top_mean    = np.mean(X_top,    axis=0)
    X_bottom_std = np.std(X_bottom, axis=0, ddof=0)
    X_top_std    = np.std(X_top,    axis=0, ddof=0)
    X_top_bottom_score = (X_top_mean-X_bottom_mean) / (X_bottom_std+X_top_std+1e-5)

    X_top_bottom_score_sorted_idx = sorted(range(len(X_top_bottom_score)), key=lambda k: abs(X_top_bottom_score[k]), reverse=True)
    X_top_bottom_score_sorted = X_top_bottom_score[X_top_bottom_score_sorted_idx]
    
    return X_top_bottom_score_sorted_idx[:n_para], X_top_bottom_score_sorted[:n_para]

# Top-bottom median filter: compare the median of gene expressions `X` of two sets of cell lines with highest and lowest `top_bottom_ratio` zscore; it is one of the dispersion filter;
# return `n_para` genes with best scores
def get_topbottom_median_filter(X, Y, n_para, top_bottom_ratio):
    Y_sorted_idx = sorted(range(len(Y)), key=lambda k: Y[k])
    
    if top_bottom_ratio>0.5:
        top_bottom_ratio = 0.5
    Y_sorted_cnt = int(Y.size*top_bottom_ratio)

    Y_label = np.zeros(Y.shape)
    Y_label[Y_sorted_idx[:Y_sorted_cnt]] = 1     # bottom
    Y_label[Y_sorted_idx[-Y_sorted_cnt:]] = 2      # top  

    X_bottom = X[Y_label==1, :]
    X_top    = X[Y_label==2, :]
    X_bottom_median = np.median(X_bottom, axis=0)
    X_top_median    = np.median(X_top,    axis=0)
    X_top_bottom_score = (X_bottom_median - X_top_median) 

    X_top_bottom_score_sorted_idx = sorted(range(len(X_top_bottom_score)), key=lambda k: abs(X_top_bottom_score[k]), reverse=True)
    X_top_bottom_score_sorted = X_top_bottom_score[X_top_bottom_score_sorted_idx]
    
    return X_top_bottom_score_sorted_idx[:n_para], X_top_bottom_score_sorted[:n_para]

# Binarize `Y` with `Y_th` and `Y_th_2`, then compare Fisher's discriminant of gene expression `X`
def get_threshold_filter_gc(X, Y, Y_th, Y_th_2=math.nan):
    if math.isnan(Y_th_2):
        Y_th_2 = Y_th
    
    Y_label = np.zeros(Y.shape)
    Y_label[Y<Y_th]   = 1    # bottom
    Y_label[Y>Y_th_2] = 2    # top
    
    X_1 = X[Y_label==1, :]
    X_2 = X[Y_label==2, :]
    X_1_mean = np.mean(X_1, axis=0)
    X_2_mean = np.mean(X_2, axis=0)
    X_1_std  = np.mean(X_1, axis=0)
    X_2_std  = np.mean(X_2, axis=0)
    X_12_score = (X_1_mean - X_2_mean) / (X_1_std+X_2_std+1e-5)
    
    X_12_score_sorted_idx = sorted(range(len(X_12_score)), key=lambda k:X_12_score[k], reverse=True)
    X_12_score_sorted = X_12_score[X_12_score_sorted_idx]
    
    # remove those lower than 0
    for i in range(len(X_12_score_sorted)):
        if X_12_score_sorted[i]<0:
            X_12_score_sorted = X_12_score_sorted[:i]
            X_12_score_sorted_idx = X_12_score_sorted_idx[:i]
            break
    
    return X_12_score_sorted_idx, X_12_score_sorted


# binarize `X` with `X_th`; Binarize `Y` with `Y_th` and `Y_th_2`; 
# return f1score
def get_double_threshold_filter_gc_f1score(X, Y, X_th, Y_th, Y_th_2=math.nan):
    X = (X > math.log(1+X_th)).astype('float')
    print(np.max(X))
    
    if math.isnan(Y_th_2):
        Y_th_2 = Y_th
    
    Y_label = np.zeros(Y.shape)
    Y_label[Y<Y_th]   = 1    # bottom
    Y_label[Y>Y_th_2] = 0    # top
    
    F1_score_list = np.zeros(X.shape[1])
    
    for i in range(X.shape[1]):
        X_sub = X[:, i]
        TP = np.sum((X_sub==1) * (Y_label==1))
        TN = np.sum((X_sub==0) * (Y_label==0))
        FP = np.sum((X_sub==1) * (Y_label==0))
        FN = np.sum((X_sub==0) * (Y_label==1))
        precision = TP / (TP + FP)
        recall = TP / (TP + FN)
        F1_score = 2 * precision * recall / (precision + recall)
        if not math.isnan(F1_score):
            F1_score_list[i] = F1_score
        
    F1_score_sorted_idx = sorted(range(len(F1_score_list)), key=lambda k:F1_score_list[k], reverse=True)
    F1_score_sorted = F1_score_list[F1_score_sorted_idx]
    
    return F1_score_sorted_idx, F1_score_sorted

# binarize `X` with `X_th`; Binarize `Y` with `Y_th` and `Y_th_2`; 
# return accuracy
def get_double_threshold_filter_gc_acc(X, Y, X_th, Y_th, Y_th_2=math.nan):
    X = (X > math.log(1+X_th)).astype('float')
    print(np.max(X))
    
    if math.isnan(Y_th_2):
        Y_th_2 = Y_th
    
    Y_label = np.zeros(Y.shape)
    Y_label[Y<Y_th]   = 1    # bottom
    Y_label[Y>Y_th_2] = 0    # top
    
    acc_list = np.zeros(X.shape[1])
    
    for i in range(X.shape[1]):
        X_sub = X[:, i]
        TP = np.sum((X_sub==1) * (Y_label==1))
        TN = np.sum((X_sub==0) * (Y_label==0))
        FP = np.sum((X_sub==1) * (Y_label==0))
        FN = np.sum((X_sub==0) * (Y_label==1))
        
        acc = (TP + TN) / (TP+TN+FP+FN)
        
        if not math.isnan(acc):
            acc_list[i] = acc
        
    acc_sorted_idx = sorted(range(len(acc_list)), key=lambda k:acc_list[k], reverse=True)
    acc_sorted = acc_list[acc_sorted_idx]
    
    return acc_sorted_idx, acc_sorted

