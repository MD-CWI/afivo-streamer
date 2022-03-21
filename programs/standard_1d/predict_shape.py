#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import interpolate
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
import joblib

p = argparse.ArgumentParser()
p.add_argument('-dataset', type=str, default='dataset.npz',
               help='Dataset filename')
p.add_argument('-N', type=int, default=10,
               help='N times low resolution')
p.add_argument('-ix', type=int, default=2,
               help='Which to reconstruct, 0 e, 1 field, 2 rhs')
args = p.parse_args()

f = np.load(args.dataset)
all_data = f['X']
ix = args.ix

# Input features
X = np.column_stack([all_data[:, 1, 0],
                     all_data[:, 1, -1]
                    ])          #, all_data[:, 0, 0]

# Target
y = all_data[:, ix, :]

X_train, X_test, y_train, y_test = \
    train_test_split(X, y, shuffle=True, random_state=1)

model = LinearRegression()

model.fit(X_train, y_train)
y_pred = model.predict(X)

fig, ax = plt.subplots(4, sharex=True)

for y, y_pred in zip(y, y_pred):
    ax[0].plot(y)
    ax[0].set_title('Data')
    ax[1].plot(y_pred)
    ax[1].set_title('Predictions')
    ax[2].plot(y_pred-y)
    ax[2].set_title('Error')

ax[3].plot(model.coef_ / model.coef_.max(), label='coef')
ax[3].plot(model.intercept_ / model.intercept_.max(), label='intercept')
ax[3].set_title('Shape function (normalized)')
ax[3].legend()

model_name = 'model_ix'+str(ix)+'_'+args.dataset[:-4]
joblib.dump(model,model_name)
print('Save model:', model_name)

n_t = all_data.shape[0]
n_shift_max = all_data.shape[2]//4
new_n_x = 3*all_data.shape[2]//4 - 10
data_new = np.zeros((n_t, all_data.shape[1], new_n_x))

ix_start = np.random.random_integers(0, n_shift_max, n_t)
for i in range(n_t):
    data_new[i,:,:] = all_data[i,:,ix_start[i]:ix_start[i]+new_n_x]

fig, ax = plt.subplots(3)
for x in data_new[:,2,:]:
    ax[0].plot(x)
print('Plot data_new rhs')

# Get a low resolution matrix
N = args.N
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
loc_orig = np.zeros(data_small.shape[0], dtype=int)

# find the location of the peak value from low resolution data
for i in range(data_small.shape[0]):
    # f = interpolate.interp1d(data_small[i, 3, :], data_small[i, 2, :],
                             # kind='quadratic', bounds_error=False)
    f = interpolate.Rbf(data_small[i, 3, :], data_small[i, 2, :],
                        bounds_error=False)
    ax[2].plot(f(data_new[i, 3, :]))
    loc_orig[i] = np.nanargmax(f(data_new[i, 3, :]))
loc_shift = all_data.shape[2]//2 - loc_orig

print('Mean shift error: ', (loc_shift-ix_start).mean())
print('Mean absolute shift error: ', np.abs(loc_shift-ix_start).mean())

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

plt.show()
