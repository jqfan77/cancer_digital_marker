import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import pickle

from read_csv_data import *
from my_filter import *
from my_load_data import * 
from my_predict import *
from my_util import *


path_root = './data/'

Y_raw, Y_drug_id, Y_drug_name, tissue_id_list, tissue_name_list, dataset_id, X_gene_symbol_list, X_oncogene_idx = my_load_data(path_root)  # 加载主要数据
cell_name_list = my_load_cell_names(path_root)

with open(path_root+'X_raw.pkl', 'rb') as f:
    X_raw = pickle.load(f)

def predict( genefile_filename, select_dataset, tissue_selected_id_list, drug_id_list ):

    # check input

    # output
    info_str = ''
    list_Y_predict = ['???' for i in range(len(drug_id_list))]
    list_Y_predict_log_ic50 = ['???' for i in range(len(drug_id_list))]
    list_Y_predict_relative = ['???' for i in range(len(drug_id_list))]
    list_color = ['#000000' for i in range(len(drug_id_list)) ]

    # 加载输入文件
    try:
        temp_rnaseq = read_csv_data_autodelimiter(genefile_filename)
        X_test_gene = []
        X_test_temp = []

        for i in range(len(temp_rnaseq)):
            if len(temp_rnaseq[i])<2:
                continue
            #基因名转数字
            if temp_rnaseq[i][0] in X_gene_symbol_list:
                X_test_gene.append( X_gene_symbol_list.index(temp_rnaseq[i][0]) )
            else:
                X_test_gene.append( 0 )
            #表达量

            X_test_temp.append( float( temp_rnaseq[i][1] ) )

        X_test = np.zeros( [1, len(X_test_temp) ] )
        X_test[0,:] = X_test_temp
        X_test = np.log(X_test + 1)
    except Exception as r:
        info_str = info_str + '[Error][file] cannot load gene expression file:\n %s\n'%(r)
        return info_str, list_Y_predict_relative, list_Y_predict, list_Y_predict_log_ic50, list_color


    # 预测时用的基因列表
    # X_test中记录的那些基因
    the_demo_gene_list = [6151, 15604, 1256, 37196, 33709, 35582, 5999, 36098, 11812, 30776,      
        5533, 5982, 15497, 30581, 33843, 4478, 2592, 12727, 5308, 8983, 
        36236, 30549, 9758, 35266, 3950, 20207, 33321, 30206, 1287, 21923, 
        34188, 33451, 30593, 30584, 835, 15881, 7854, 2035, 11894, 30589, 
        7444, 19629, 874, 24443, 836, 2053, 30000, 23749, 32956, 22450, 
        2012, 10987, 8656, 30201, 34808, 24234, 873, 937, 2989, 19294, 
        30250, 30503, 36293, 3037, 33271, 15385, 9317, 28554, 23322, ]
    the_gene_list_for_test_dict = { 
        1005: the_demo_gene_list, 
        1073: the_demo_gene_list,
        1080: the_demo_gene_list,
        1512: the_demo_gene_list
    }
    C_neighbor_k_dict = {
        1005: 1, 
        1073: 1,
        1080: 1,
        1512: 1
    }


    # 准备开始
    

    for iii in range(len(drug_id_list)):

        drug_id = drug_id_list[iii]

        if not drug_id in the_gene_list_for_test_dict:
            return '[Error] Unknown drug_id %d'%(drug_id)
        X_idx = the_gene_list_for_test_dict[drug_id]
        C_neighbor_k = C_neighbor_k_dict[drug_id]
        drug_idx = np.argwhere(Y_drug_id==drug_id)[0][0]

        # 最终决定X_idx
        X_idx = np.array(X_idx)
        f_idx = get_repeating_idx_filter(X_idx, np.array(X_test_gene))
        if f_idx.shape[0]<10:
            info_str = info_str + '[Error][drug%d] Not enough gene expression for test. \nThere are only %d genes available for the test\n' % (drug_id, f_idx.shape[0])
            continue
        X_idx = X_idx [f_idx]


        # 获得test的提取表
        X_test_idx = np.zeros(X_idx.shape, np.int)
        for j in range(X_idx.shape[0]):
            X_test_idx[j] = X_test_gene.index(X_idx[j])


        # 开始提取数据进行预测
        X = X_raw[:, X_idx]
        X = np.log(1+X)
        Y = Y_raw[drug_idx, :]
        # 按照数据库、组织，选择部分细胞
        X_cell_idx = my_data_selection_2(dataset_id, select_dataset, tissue_id_list, tissue_selected_id_list, False)
        # 删除没药敏的
        X_cell_idx = X_cell_idx[~np.isnan(Y[X_cell_idx])]
        # 应用
        X = X[X_cell_idx, :]
        Y = Y[X_cell_idx]
        # 模型
        model_knn = KNeighborsRegressor(n_neighbors = C_neighbor_k, weights = 'uniform', algorithm = 'auto')
        model_knn.fit(X,Y)

        
        # 处理test数据，获得zscore
        X_test_copy = X_test.copy()
        X_test_copy = X_test_copy[:, X_test_idx]
        Y_predict = model_knn.predict(X_test_copy)[0]

        # 转化成log10_IC50
        the_zscore_to_log_ic50_coeff_1_dict = {
            1005:1.9455148866575,
            1080:2.07199265547472,
            1512:0.940400101814217,
            1073:1.71149424488917,
        }
        the_zscore_to_log_ic50_coeff_2_dict = {
            1005:3.62032026591016,
            1080:-2.55823666134551,
            1512:5.3947264425877,
            1073:4.6811504770386,
        }
        coeff_1 = the_zscore_to_log_ic50_coeff_1_dict[drug_id]
        coeff_2 = the_zscore_to_log_ic50_coeff_2_dict[drug_id]
        Y_predict_log_ic50 = ( Y_predict*coeff_1 + coeff_2) / math.log(10)

        # 转化成相对zscore
        Y_mean = np.mean(Y)
        Y_std = np.std(Y)
        Y_predict_relative = ( Y_predict - Y_mean ) / Y_std
        # 并获得对应的颜色
        the_color = '#000000'
        if Y_predict_relative<-1:
            the_color = '#%02X%02X%02X' % (34,177,76)
        elif Y_predict_relative>1:
            the_color = '#%02X%02X%02X' % (237,28,36)
        

        #
        list_Y_predict[iii] = '%.3f' % (Y_predict)
        list_Y_predict_log_ic50[iii] = '%.3f' % (Y_predict_log_ic50)
        list_Y_predict_relative[iii] = '%.3f' % (Y_predict_relative)
        list_color[iii] = the_color

    return info_str, list_Y_predict_relative, list_Y_predict, list_Y_predict_log_ic50, list_color
