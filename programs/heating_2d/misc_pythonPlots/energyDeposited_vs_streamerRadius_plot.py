
import numpy as np
import matplotlib.pyplot as plt

#Loading the r_vs_z.txt file

rData = np.loadtxt("r_vs_z.txt")
#Shifting the origing to the top of the domain
rData[:,0] = 16e-3 - rData[:,0]
#Loading the axial energy distribution file

eData = np.loadtxt("line_energy_density_0054.curve")

e_interp = np.interp(rData[:,0], eData[:,0], eData[:,1])

fig, ax = plt.subplots()
ax.plot(rData[1:-1,1], e_interp[1:-1],'r.-',label='J.E')
ax.plot(rData[1:-1,1], rData[1:-1,1]**(-4)*5e-10,'b.-', label='K/$r^4$')
ax.legend()
plt.show()
