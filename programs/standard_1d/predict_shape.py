#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

p = argparse.ArgumentParser()
p.add_argument('-dataset', type=str, default='dataset.npz',
               help='Dataset filename')
args = p.parse_args()

f = np.load(args.dataset)
all_data = f['X']

# Input features
X = np.column_stack([all_data[:, 1, 0],
                     all_data[:, 1, -1],
                     all_data[:, 0, 0]])

# Target
y = all_data[:, 1, :]

X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True)

# pca = PCA(n_components=5).fit(y)
# y = pca.transform(y)

model = LinearRegression()
# model = MultiOutputRegressor(GradientBoostingRegressor())
# model = make_pipeline(StandardScaler(), MLPRegressor())

model.fit(X_train, y_train)
y_pred = model.predict(X_test)
# y_pred = cross_val_predict(model, X, y)

# y_pred = pca.inverse_transform(y_pred)
# y = pca.inverse_transform(y)

fig, ax = plt.subplots(4, sharex=True)

model.fit(X, y)

for y, y_pred in zip(y_test, y_pred):
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

plt.show()
