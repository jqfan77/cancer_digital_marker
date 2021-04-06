import numpy as np
import csv
import os
from read_csv_data import *

## load some basic information
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

    ### load Y
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
    # X_gene_symbol； # X_oncogene_symbol 
    X_gene_symbol_list = []
    for i in range(len(X_gene_symbol)):
        X_gene_symbol_list.append(X_gene_symbol[i][2])
    X_oncogene_idx = []
    for i in range(len(X_oncogene_symbol)):
        if X_oncogene_symbol[i][0] in X_gene_symbol_list:
            X_oncogene_idx.append( X_gene_symbol_list.index(X_oncogene_symbol[i][0]) )

    return Y_raw, Y_drug_id, Y_drug_name, tissue_id_list, tissue_name_list, dataset_id, X_gene_symbol_list, X_oncogene_idx

## load matrix X
def my_load_data_X(path):
    X_raw = read_csv_data(path, firstline=False)
    X_raw = [[float(x) for x in R]for R in X_raw]
    X_raw = np.array(X_raw)
    return X_raw

###########################################################
## load data from TCGA
# select data of the specific drug
def my_select_drug_name_tcga(path_drug_sens,drug_selected_name='ALL'):
    y_1=read_csv_data(path_drug_sens, row=-1,firstline=False)
    L=len(y_1)
    uuid=[]
    drug_sens=[]
    if drug_selected_name == 'ALL':
        uuid=np.array(y_1)[:,0]
        drug_sens=np.array(y_1)[:,2].astype(np.float)
    else:
        for i in range(L):
            if y_1[i][1]==drug_selected_name and y_1[i][0] not in uuid:
                uuid.append( y_1[i][0])
                drug_sens.append(float(y_1[i][2]))
        uuid=np.array(uuid)
        drug_sens=np.array(drug_sens)
    return uuid,drug_sens

# match the case uuid and filenames
def my_match_uuid_filename_tcga(path_drug_match,u_1,y_1):
    match_file=read_csv_data(path_drug_match, row=-1,firstline=True)
    L_1=len(u_1)
    L_2=len(match_file)
    f=[]
    flag_not_null=[]
    for i in range(L_1):
        flag=0
        for j in range(L_2):
            if u_1[i]==match_file[j][0]:
                f.append(match_file[j][1])
#                 print(f[i]+'     over')
                flag=1
                break
        if flag==0:
            f.append('')
        if len(f[i])!=0:
            flag_not_null.append(i)
    f=np.array(f)        
    f=f[flag_not_null]
    u_1=u_1[flag_not_null]
    y_1=y_1[flag_not_null]
        
    return f,u_1,y_1

# select RNA-seq data
def my_select_X_tcga(path_root,X_raw,f_1,u_1,y_1):

    dirname = path_root+'/FPKM/FPKM.txt/'
    filename_list = os.listdir(dirname)

    filename_cnt = len(filename_list)
    L=len(f_1)
    X=np.zeros([L,np.size(X_raw,1)],dtype=float)
    flag_match=[]
    
    for i in range(L):
        flag=0
        for j in range(filename_cnt):
            if f_1[i]==filename_list[j]:
                X[i]=X_raw[j]
                flag=1
                break
        if flag==0:
            print('FKPM file',i,'mismatch','deleted')
        else:
            flag_match.append(i)
   
    X=X[flag_match]
    f_1=f_1[flag_match]
    u_1=u_1[flag_match]
    y_1=y_1[flag_match]  
    return X,f_1,u_1,y_1

# delete genes without expression
def my_delete_NaN_col(X1):
    X=X1.copy()
    num_col=np.size(X,1)
    flag_not_NaN_col=[]
    for i in range(num_col):
        flag=0
        for j in range(np.size(X[:,i])):
            if X[j,i]<0:
                flag=1
                break
        if flag==0:
            flag_not_NaN_col.append(i)
    return X[:,flag_not_NaN_col],flag_not_NaN_col

# load TCGA data
def my_load_data_tcga(tissue_selected,drug_name):
    
    path_root = '../../TCGA-'+tissue_selected+'/'
    path_X = path_root + 'FPKM/FPKM_TCGA_'+tissue_selected+'_GDSCfmt.csv'
    path_y_sens = path_root + 'drug_'+tissue_selected+'/drug_'+tissue_selected+'_sens.csv'
    path_y_match = path_root + 'drug_'+tissue_selected+'/drug_'+tissue_selected+'_match.csv'
    
    X_raw = my_load_data_X(path_X)
    # select the specific drug, return uuid list and drug sensitivity list y
    u_1,y_1=my_select_drug_name_tcga(path_y_sens,drug_name)
    # fpkm filename list，uuid list and y
    f_1,u_1,y_1=my_match_uuid_filename_tcga(path_y_match,u_1,y_1)
    # obtain fpkm array
    X,f_1,u_1,y_1=my_select_X_tcga(path_root,X_raw,f_1,u_1,y_1)
    
    return X,f_1,u_1,y_1
