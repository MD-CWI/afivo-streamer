#!/usr/bin/env/ python



import argparse
import os
import numpy as np
import matplotlib.pyplot as plt



pr = argparse.ArgumentParser(
        description='''Plot a variable at a point from line data for each time step.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: python3 plot_point_time.py -dbname database -dbloc foldername -point pt -timeStep dt''')
# Improvement: we can give multiple points and then we get multiple graphs in the same plot
pr.add_argument('-dbname', type=str, nargs=1, required=True,
                help='Database name (eg line_varname(excluding the .curve))')
pr.add_argument('-dbloc', type=str, nargs=1, required=True,
                help='Database location folder name')
pr.add_argument('-point', type=float, nargs=1, required=True,
                help='Point at which the data is taken')
pr.add_argument('-timeStep', type=float, nargs=1, required=True,
                help='Time step so as to calculate the time')
#cmap = plt.cm.get_cmap('plasma')
#cmap = plt.cm.get_cmap('hot')
#linecolor = cmap(np.linspace(0,1,100)) #Limiting here the number of lines/linecolors to 100 as we do not want more than 100 plots in a graph
#line_iter = 0
args = pr.parse_args()
fig, ax = plt.subplots()
for fname in os.listdir(args.dbloc[0]):
    if fname.startswith(args.dbname[0]):
        time = int(fname.split('_')[2].split('.')[0])*args.timeStep[0]
        t = np.loadtxt(args.dbloc[0]+fname, delimiter=' ')
        y = np.interp(args.point[0], t[:,0], t[:,1])
        #ax.plot(t[:,0], t[:,1], c=linecolor[line_iter])
        ax.plot(time, y, 'bo')



plt.show()
