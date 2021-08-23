# Last updated on Aug 21th, 2021,  fjq19@mails.tsinghua.edu.cn
import numpy as np

# normalization of the matrix X
def my_normalize(X):
    if X.ndim==1:
        X_mean = np.nanmean(X)
        X_std = np.nanstd(X)
        X_normalized = (X - X_mean) / X_std
        return X_normalized
    elif X.ndim==2:
        X_mean = np.nanmean(X, axis=0)
        X_std = np.nanstd(X, axis=0)
        X_std[X_std==0] = 1
        X_normalized = (X - X_mean) / X_std
        return X_normalized
    else:
        X_mean = np.nanmean(X)
        X_std = np.nanstd(X)
        X_normalized = (X - X_mean) / X_std
        return X_normalized
    
# TPM normalization
def fpkm_to_tpm(x):
    x_sum = np.sum(x, axis=1) / 1e6
    for i in range(x.shape[0]):
        x[i,:] = x[i,:] / x_sum[i]
    return x