#!/usr/bin/env python

import numpy as np
import collections
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

class Plotdata():

    def __init__(self):
        self.d = collections.OrderedDict()

    def add_file(self, fname, label):
        self.d[label] = np.genfromtxt(fname, names=True)

    def plotz(self):
        plt.clf()
        for key in self.d.keys():
            t = self.d[key]
            plt.plot(t[self.xvar], t['z_Ez'])
        plt.legend(self.d.keys())
        plt.ylabel('z (m)')
        plt.xlabel(self.xvar)
        plt.show()

    def plot(self, xname, yname, filter_width=0):
        plt.clf()
        for key in self.d.keys():
            t = self.d[key]
            y = t[yname]
            if filter_width > 0:
                y = savgol_filter(dy, 2 * filter_width + 1, 2)
            plt.plot(t[xname], y)
        plt.legend(self.d.keys())
        plt.ylabel(yname)
        plt.xlabel(xname)
        plt.show()

    def plot_deriv(self, xname, yname, filter_width=0):
        plt.clf()
        for key in self.d.keys():
            t = self.d[key]
            dy = np.diff(t[yname])/np.diff(t[xname])
            if filter_width > 0:
                dy = savgol_filter(dy, 2 * filter_width + 1, 2)
            plt.plot(t[xname][:-1], dy)
        plt.legend(self.d.keys())
        plt.ylabel(yname)
        plt.xlabel(xname)
        plt.show()

if __name__ == '__main__':
    p = Plotdata()

    p.add_file("output/guiding_cyl.txt", "test")
    p.add_file("output/guiding_cyl_2.txt", "test2")
    # p.add_file("output/guiding_cyl_adx6.txt", "a6")
    p.add_file("output/guiding_cyl_adx8.txt", "a8")
    # p.add_file("output/guiding_cyl_adx10.txt", "a10")
    # p.add_file("output/guiding_cyl_adx12.txt", "a12")
    # p.add_file("output/guiding_cyl_adx14.txt", "a14")
    # p.add_file("output/guiding_cyl_adx16.txt", "a16")
    # p.plot('time', 'z_Ez')
    # p.plot('time', 'r_Er')
    # p.plot_deriv('time', 'Er', 2)
    p.plot_deriv('time', 'z_Ez')
