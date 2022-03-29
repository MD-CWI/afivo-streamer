#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 6 18:22:07 2021

@author: Baohong.Guo
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from my_module import generate_y_field

font = {
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Times New Roman",
    "font.size": 14
    }
mpl.rcParams.update(font)

fig, axes = plt.subplots(1, 1, figsize=(6, 4))
plt.subplots_adjust(wspace=0.15, hspace=0.04)

L = 0.01
y_field = generate_y_field(profile='exponential', initial_field=2e6, goal_field=1e6, 
                     decay_distance=5e-3, z_constant=5e-3, z_fall=95e-3)
y = y_field[:,0]
field = y_field[:,1]
voltage = y_field[:,1]*L

ax = axes
ax.plot(y/1e-3, field/1e5, label=r'$E_\mathrm{bg}$')

ax.set_xlabel(r'$y$ (mm)')
ax.set_ylabel(r'$E_\mathrm{bg}$ (kV/cm)')
#ax.set_xlim(50, 0)
#ax.set_ylim(2.5, 13.5)
ax.minorticks_on()
ax.tick_params(top=True, bottom=True, left=True, right=True)
ax.tick_params(axis='both', which='both', direction='in')
ax.grid(axis='both', which='major', linestyle='-.', linewidth=0.5, alpha=0.7)
ax.legend(fontsize=12, ncol=1, borderaxespad=0.3, columnspacing=0.5, labelspacing=0.3, handletextpad=0.5, frameon=False)

# =============================================================================
# np.savetxt('y_field_voltage.txt', np.column_stack([y, field, voltage]), 
#            header='y  field  voltage', comments='', delimiter='  ', fmt='%.8e')  
# print('Saved y_field_voltage.txt')
# =============================================================================

#plt.savefig('generate_y_field_voltage.png', dpi=600, bbox_inches='tight', pad_inches=0.05)
#plt.savefig('generate_y_field_voltage.pdf', bbox_inches='tight', pad_inches=0.05)
plt.show()

