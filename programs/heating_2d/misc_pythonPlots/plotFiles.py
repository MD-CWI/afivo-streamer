#!/usr/bin/env/ python
# Plots the .curve files for all time steps of a given variable.


import argparse
import os
import numpy as np
import matplotlib.pyplot as plt



pr = argparse.ArgumentParser(
        description='''Plot the set of files provided''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: python3 plotFiles.py file1 file2 ...''')

pr.add_argument("data_files", type=str, nargs='+', help="File(s) to plot")
args = pr.parse_args()

files = [np.loadtxt(f) for f in args.data_files]

# Add functionality for different graph colors/shapes
cmap = plt.cm.get_cmap('gist_rainbow')

#fig,ax = plt.subplots()

for i,f in enumerate(files):
    c = cmap(float(i)/ len(files))
    plt.plot(f[:,0], f[:,1], color = c, label = args.data_files[i])
plt.legend()
plt.show()
