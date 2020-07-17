import numpy as np
import matplotlib.pyplot as plt

filename = 'streamer_2d_log.txt'

#Loading the file contents as a 2 dimensional numpy array

baohong = np.loadtxt(filename, skiprows= 1) #here we skip the first row as it is just text


#Check what we need to plot. Make sure you have consistent units for the quantities
# Here, I am trying to plot L-vt

t = baohong[:,1]
#v = baohong[:,3]
L = baohong[:,8]

plt.plot(t/1e-8, L*1e2 - 0.03*t, 'bo-')
plt.xlabel('t[ns]')
plt.ylabel('L(t) - vt[cm]')
plt.show()

