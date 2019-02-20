#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 09:04:50 2019

Testing correlation code in hail.

@author: nbaya
"""

#import hail as hl
from hail.linalg import BlockMatrix as bmat
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


X_np = np.random.normal(size=[10,30])
X_np[:,:int(X_np.shape[1]*0.8)] = np.random.normal(loc=2,scale=0.000000000000001,size=X_np[:,:int(X_np.shape[1]*0.8)].shape)
X_np[:,:int(X_np.shape[1]*0.2)] = np.random.normal(loc=-2,scale=0.00000000000001,size=X_np[:,:int(X_np.shape[1]*0.2)].shape)
X_np[:,:int(X_np.shape[1]*0.1)] = np.random.normal(loc=0,scale=0.0000000000000000001,size=X_np[:,:int(X_np.shape[1]*0.1)].shape)
#X_np[:int(X_np.shape[0]*0.5),:] = np.random.normal(loc=0,scale=0.001,size=X_np[:int(X_np.shape[0]*0.5),:].shape)
#X_np[:int(X_np.shape[0]*0.4),:] = np.random.normal(loc=2,scale=0.001,size=X_np[:int(X_np.shape[0]*0.4),:].shape) #np.random.random(size=X_np[:int(X_np.shape[0]*0.4),:].shape)
#X_np[:int(X_np.shape[0]*0.2),:] = np.random.exponential(size=X_np[:int(X_np.shape[0]*0.2),:].shape)

R_np = np.corrcoef(X_np.T)

X_hl = bmat.from_numpy(X_np)
X_hl.shape
CX = X_hl-X_hl.sum(axis=0)/X_hl.shape[0] #mean-centered matrix
cov = (1/X_hl.shape[0])*CX.T@CX
D = (cov**(-1/2)).sparsify_band(lower=0,upper=0)
X_s = CX@D
R = (1/X_hl.shape[0])*X_s.T@X_s
R.shape

R_hl = R.to_numpy()

print('normed difference: '+str(np.linalg.norm(R_hl-R_np)))

fig, ax = plt.subplots()
sns.heatmap(R_hl,cmap='RdYlGn',ax=ax)
fig.set_size_inches(12,10)
fig=plt.gcf()