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
from my_module import generate_t_field

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
t_field = generate_t_field(tStart=5.0e-9, tRise=5.0e-9, tPulse=10.0e-9, tFall=5.0e-9, tPost=10.0e-9, \
                     field_start=1.0e6, field_pulse=2.0e6, field_post=0.0e6)
t = t_field[:,0]
field = t_field[:,1]
voltage = t_field[:,1]*L

ax = axes
ax.plot(t/1e-9, field/1e5, label=r'$E_\mathrm{bg}$')

ax.set_xlabel(r'time (ns)')
ax.set_ylabel(r'$E_\mathrm{bg}$ (kV/cm)')
#ax.set_xlim(50, 0)
#ax.set_ylim(2.5, 13.5)
ax.minorticks_on()
ax.tick_params(top=True, bottom=True, left=True, right=True)
ax.tick_params(axis='both', which='both', direction='in')
ax.grid(axis='both', which='major', linestyle='-.', linewidth=0.5, alpha=0.7)
ax.legend(fontsize=12, ncol=1, borderaxespad=0.3, columnspacing=0.5, labelspacing=0.3, handletextpad=0.5, frameon=False)

# =============================================================================
# np.savetxt('time_field_voltage.txt', np.column_stack([t, field, voltage]), 
#            header='time  field  voltage', comments='', delimiter='  ', fmt='%.8e')  
# print('Saved y_field_voltage.txt')
# =============================================================================

#plt.savefig('generate_time_field_voltage.png', dpi=600, bbox_inches='tight', pad_inches=0.05)
#plt.savefig('generate_time_field_voltage.pdf', bbox_inches='tight', pad_inches=0.05)
plt.show()

