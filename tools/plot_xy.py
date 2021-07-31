#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 10:18:17 2021

@author: Baohong.Guo
"""

import numpy as np
import pandas as pd
import argparse
import sys
import os
import matplotlib.pyplot as plt

# =============================================================================
# p = argparse.ArgumentParser(
#     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# 
# p.add_argument("log_file", type=str, nargs= '+', help="Input log file(s)")
# p.add_argument("x", type=str, default='time', help="Name of x variable")
# p.add_argument("y", type=str, default='max(E)', help="Name of y variable")
# args = p.parse_args()
# 
# logs = [pd.read_csv(f, delim_whitespace=True) for f in args.log_file]
# numbered_files = [f'{i}: {f}' for i, f in enumerate(args.log_file)]
# 
# fig, axes = plt.subplots(1, 1, constrained_layout=True)
# fig.suptitle('\n'.join(numbered_files))
# 
# for i, log in enumerate(logs):
#     axes.plot(log[args.x], log[args.y], label=f'{args.y}-{i}')
#     
# plt.xlabel(args.x)
# plt.ylabel(args.y)
# plt.legend()
# plt.show()
# =============================================================================

x, y = input('Enter names of variable x y:').split()

filepath = r'C:\Users\Baohong.Guo\Desktop\3'
filename = os.listdir(filepath)

logs = [pd.read_csv(f, delim_whitespace=True) for f in filename]
numbered_files = [f'{i}: {f}' for i, f in enumerate(filename)]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    axes.plot(log[x], log[y], label=f'{y}-{i}')

plt.xlabel(x)
plt.ylabel(y)
plt.legend()
plt.show()

