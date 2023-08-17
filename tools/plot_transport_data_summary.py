#!/usr/bin/env python3

# Author: Jannis Teunissen
# Purpose: plot the _summary.txt file from a simulation

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument("summary_file", type=str, help="File <simulation>_summary.txt")
p.add_argument("-SI_field", action='store_true',
               help="Use electric field in V/m rather than Td")
args = p.parse_args()

if args.SI_field:
    tdata = pd.read_csv(args.summary_file, delim_whitespace=True, index_col=1)
    tdata.drop(columns='E/N[Td]', inplace=True)
else:
    tdata = pd.read_csv(args.summary_file, delim_whitespace=True, index_col=0)
    tdata.drop(columns='E[V/m]', inplace=True)

fig, ax = plt.subplots(figsize=(10, 10))
tdata.plot(subplots=True, layout=(-1, 2), ax=ax, sharex=True)

plt.show()
