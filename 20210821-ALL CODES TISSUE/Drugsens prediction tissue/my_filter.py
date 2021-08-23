# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn
import numpy as np
from scipy.stats import mannwhitneyu

# mannwhitneyu-test
def get_utest_filter(X,Y,th_y,th = 0.1):
    Filter_idx = []
    Filter_pvalue = []
    for i in range(X.shape[1]):
        X_flag = X[:,i]
        x_p = mannwhitneyu(X_flag[Y>th_y],X_flag[Y<=th_y])[1]
        if x_p<th:
            Filter_idx.append(i)
            Filter_pvalue.append(x_p)
    Filter_pvalue=np.array(Filter_pvalue)
    return Filter_idx, Filter_pvalue

# select genes in idx_list_2 from idx_list_1
def get_repeating_idx_filter(idx_list_1, idx_list_2):
    Filter_idx = []
    for idx_2 in idx_list_2:
        if idx_2 in idx_list_1:
            Filter_idx.append(np.argwhere(idx_list_1==idx_2)[0][0])
    return np.array(Filter_idx)