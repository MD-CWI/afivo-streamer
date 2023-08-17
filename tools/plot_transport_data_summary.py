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

tdata = pd.read_csv(args.summary_file, delim_whitespace=True)

if args.SI_field:
    tdata.set_index('E[V/m]', inplace=True)
    tdata.drop(columns='E/N[Td]', inplace=True)
else:
    tdata.set_index('E/N[Td]', inplace=True)
    tdata.drop(columns='E[V/m]', inplace=True)

tdata.plot(subplots=True, layout=(-1, 2), sharex=True, figsize=(10, 10))

plt.show()
