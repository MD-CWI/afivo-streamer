#!/usr/bin/env/ python



import argparse
import os
import numpy as np
import matplotlib.pyplot as plt



pr = argparse.ArgumentParser(
        description='''Plot a set of line data for each time step.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: python3 plot_line_time.py -dbname database -dbloc foldername''')

pr.add_argument('-dbname', type=str, nargs=1, required=True,
                help='Database name (eg line_varname(excluding the .curve))')
pr.add_argument('-dbloc', type=str, nargs=1, required=True,
                help='Database location folder name)')

args = pr.parse_args()
fig, ax = plt.subplots()
for fname in os.listdir(args.dbloc[0]):
    if fname.startswith(args.dbname[0]):
        t = np.loadtxt(args.dbloc[0]+fname, delimiter=' ')
        ax.plot(t[:,0], t[:,1])



plt.show()
