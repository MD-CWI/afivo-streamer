#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 12:27:18 2022

@author: xiaoran
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import interpolate
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import joblib

p = argparse.ArgumentParser()
p.add_argument('-dataset', type=str, default='dataset.npz',
               help='Dataset filename')
p.add_argument('-model', type=str, default='model_dataset',
               help='Model filename')
p.add_argument('-ix', type=int, default=2,
               help='Which to reconstruct, 0 e, 1 field, 2 rhs')
args = p.parse_args()

f = np.load(args.dataset)
all_data = f['X']
ix = args.ix

model = joblib.load(args.model)

n_t = all_data.shape[0]
n_shift_max = 300
new_n_x = all_data.shape[2]-n_shift_max-50
data_new = np.zeros((n_t, all_data.shape[1], new_n_x))

ix_start = np.random.random_integers(0, n_shift_max, n_t)
for i in range(n_t):
    data_new[i,:,:] = all_data[i,:,ix_start[i]:ix_start[i]+new_n_x]

fig, ax = plt.subplots(2)
for x in data_new[:,2,:]:
    ax[0].plot(x)
print('Plot data_new rhs')

# Get a low resolution matrix
N = 10
data_small = np.zeros((data_new.shape[0], data_new.shape[1], int(data_new.shape[2]/N)))
for i in range(int(data_new.shape[2]/N)):
    data_small[:,:,i] = np.average(data_new[:,:,N*i:N*i+N],axis=2)

for x in data_small[:,2,:]:
    ax[1].plot(x)
print('Plot data_small rhs')

'''
fig, ax = plt.subplots(2)
for a, x in zip (data_new[:,3,:],data_new[:,2,:]):
    ax[0].plot(a,x)
for a, x in zip (data_small[:,3,:],data_small[:,2,:]):
    ax[1].plot(a,x)
'''


X_small = np.column_stack([data_small[:, 1, 0],
                         data_small[:, 1, -1]                      
                        ]) #  , data_small[:, 0, 0]
y_small = model.predict(X_small)
#y_small = pca.inverse_transform(model.predict(X_small))

'''
plt.figure()
for y in y_small:
    plt.plot(y)
print('done')
'''

loc_pred = np.argmax(y_small, axis=1)
loc_orig = np.zeros(data_small.shape[0])
# find the location of the peak value from low resolution data
for i in range(data_small.shape[0]):
    #f = interpolate.interp1d(data_small[i, 3, :],data_small[i, 2, :], kind='quadratic',bounds_error=False)
    f = interpolate.Rbf(data_small[i, 3, :],data_small[i, 2, :],bounds_error=False)
    loc_orig[i] = np.nanargmax(f(data_new[i, 3, :]))
loc_shift = all_data.shape[2]/2 - loc_orig
loc_shift = loc_shift.astype(int)


fig, ax = plt.subplots(4)
for y, loc, loc2, y2 in zip(y_small, loc_shift, ix_start, data_new[:,ix,:]):#, data_new[:,3,:], ix_start
    ax[0].plot(y[loc:loc+new_n_x])
    ax[0].set_title('Predictions')
    ax[1].plot(y2)
    ax[1].set_title('Data')
    ax[2].plot(y2-y[loc:loc+new_n_x])
    ax[2].set_title('Error_with_predict_loction')
    ax[3].plot(y2-y[loc2:loc2+new_n_x])
    ax[3].set_title('Error_with_original_loction')
print('done')
