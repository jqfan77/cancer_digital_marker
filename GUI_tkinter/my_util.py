# Last updated on 2021 Mar 30, 657947836@qq.com


import math
import numpy as np
import matplotlib.pyplot as plt

# to ensure mean=0 and std=1
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
