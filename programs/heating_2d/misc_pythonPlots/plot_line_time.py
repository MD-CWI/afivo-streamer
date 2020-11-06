#!/usr/bin/env/ python
# Plots the .curve files for all time steps of a given variable.


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
#cmap = plt.cm.get_cmap('plasma')
cmap = plt.cm.get_cmap('hot')
linecolor = cmap(np.linspace(0,1,100)) #Limiting here the number of lines/linecolors to 100 as we do not want more than 100 plots in a graph
line_iter = 0
args = pr.parse_args()
fig, ax = plt.subplots()
print(args.dbname)
for fname in os.listdir(args.dbloc[0]):
    if fname.startswith(args.dbname[0]):
        if line_iter%10 == 0 :
        #if (fname.endswith('0000.curve') or fname.endswith('0010.curve') or fname.endswith('0020.curve') or fname.endswith('0030.curve') or fname.endswith('0040.curve') or fname.endswith('0055.curve')):
            t = np.loadtxt(args.dbloc[0]+fname, delimiter=' ')
            ax.plot(t[:,0], t[:,1])
            #ax.plot(t[:,0], t[:,1], c=linecolor[line_iter])
            print(line_iter)
        line_iter = line_iter + 1



plt.show()
