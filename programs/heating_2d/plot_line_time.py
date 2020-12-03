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
#cmap = plt.cm.get_cmap('plasma')
cmap = plt.cm.get_cmap('hot')
linecolor = cmap(np.linspace(0,1,100)) #Limiting here the number of lines/linecolors to 100 as we do not want more than 100 plots in a graph
line_iter = 0
args = pr.parse_args()
fig, ax = plt.subplots()
for fname in os.listdir(args.dbloc[0]):
    if fname.startswith(args.dbname[0]):
        if line_iter%5 == 0 :
            t = np.loadtxt(args.dbloc[0]+fname, delimiter=' ')
            ax.plot(t[:,0], t[:,1], c=linecolor[line_iter])
            print(line_iter)
        line_iter = line_iter + 1



plt.show()
