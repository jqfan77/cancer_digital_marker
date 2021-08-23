# Load data. return `np.array` for number matrix, return `list` for strings 
# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn

import numpy as np
import csv
from read_csv_data import *

# Load data
def my_load_data(path_root):
    path_Y='rnaseq_zscore_all.csv'  ##sample-drugsens

    path_gene_symbol='rnaseq_gene_symbol.csv'  ##gene name
    path_oncogene_symbol = 'rnaseq_oncogene_symbol.csv'
    path_drug_name='drug_name.csv'  ##drug name
    path_drug_id = 'drug_id.csv'    ##drug id

    path_tissue_id='cell_dataset_tissueid.csv' ##cell-tissue
    path_tissue_name='tissue_name.csv'  ##tissue name

    path_dataset_id='rnaseq_fpkm_model_dataset.csv' ## dataset 1_ccle，2_gdsc

    path_cell_name='rnaseq_fpkm_model_name.csv'
    
    
    #################################
    ### load Y, drug zscore
    Y_raw = read_csv_data(path_root+path_Y)
    Y_raw = [[float(x) for x in R]for R in Y_raw]
    Y_raw = np.array(Y_raw)
    
    ### load drug id
    Y_drug_id = read_csv_data(path_root+path_drug_id)
    Y_drug_id = [int(x[0]) for x in Y_drug_id]
    Y_drug_id = np.array(Y_drug_id)
    ### load drug names
    Y_drug_name = read_csv_data(path_root+path_drug_name)
    Y_drug_name = [x[0] for x in Y_drug_name]
    
    ### load tissue names
    tissue_name_list=read_csv_data(path_root+path_tissue_name)
    tissue_name_list=np.array(tissue_name_list)
    tissue_name_list=tissue_name_list.T[1].tolist()
    ### load tissue id of cells
    tissue_id_list=read_csv_data(path_root+path_tissue_id,firstline=False)
    tissue_id_list=np.array(tissue_id_list)
    tissue_id_list=tissue_id_list[:,2].T
    tissue_id_list = np.array([int(x) for x in tissue_id_list])

    ### load dataset id of cells
    dataset_id=read_csv_data(path_root+path_dataset_id)
    dataset_id=np.array([[int(x) for x in R]for R in dataset_id])
    dataset_id=dataset_id.T[0]

    ### load gene symbol of rnaseq
    X_gene_symbol = read_csv_data(path_root+path_gene_symbol,firstline=False)
    X_oncogene_symbol = read_csv_data(path_root + path_oncogene_symbol)
    # transform `X_oncogene_symbol` to index
    X_gene_symbol_list = []
    for i in range(len(X_gene_symbol)):
        X_gene_symbol_list.append(X_gene_symbol[i][2])
    X_oncogene_idx = []
    for i in range(len(X_oncogene_symbol)):
        if X_oncogene_symbol[i][0] in X_gene_symbol_list:
            X_oncogene_idx.append( X_gene_symbol_list.index(X_oncogene_symbol[i][0]) )

    return Y_raw, Y_drug_id, Y_drug_name, tissue_id_list, tissue_name_list, dataset_id, X_gene_symbol_list, X_oncogene_idx

################################################################################

# load X, gene expression 
def my_load_data_X(path_root):
    path_X='rnaseq_fpkm_trans.csv'  
    X_raw = read_csv_data(path_root+path_X, firstline=False)
    X_raw = [[float(x) for x in R]for R in X_raw]
    X_raw = np.array(X_raw)
    
    return X_raw

# load X, and TPM normalize
def my_load_data_X_tpm(path_root):
    path_X='rnaseq_fpkm_trans.csv'  
    X_raw = read_csv_data(path_root+path_X, firstline=False)
    X_raw = [[float(x) for x in R]for R in X_raw]
    X_raw = np.array(X_raw)
    
    X_raw_sum = np.sum(X_raw, axis=1) / 1e6
    for i in range(X_raw.shape[0]):
        X_raw[i,:] = X_raw[i,:] / X_raw_sum[i]
    return X_raw

################################################################################

# load names of cell lines
def my_load_cell_names(path_root):
    path_tissue_id='cell_dataset_tissueid.csv' ##cell-tissue
    tissue_id_list=read_csv_data(path_root+path_tissue_id,firstline=False)
    tissue_id_list=np.array(tissue_id_list)
    cell_name_list=tissue_id_list[:,0].T
    
    return cell_name_list.tolist()

# load ENSG## for genes
def my_load_ENSG(path_root):
    path_ENSG = 'rnaseq_gene_symbol_ENSG.csv' ##cell-tissue
    X_gene_ENSG_list = read_csv_data(path_root+path_ENSG,firstline=True)
    
    X_gene_ENSG_list = np.array(X_gene_ENSG_list)
    X_gene_ENSG_list = X_gene_ENSG_list[:,3].T
    
    return X_gene_ENSG_list.tolist()
