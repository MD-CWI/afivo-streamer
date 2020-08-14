import numpy as np
import matplotlib.pyplot as plt




#Electric field amplitude (Volts)
field_amplitude = -10.0e+5 

#Pulse start time
tStart = 0.0e-9

#Pulse duration (seconds)
tPulse = 30.0e-9


#Rise time (seconds)
tRise = 5.0e-9

#Fall time (seconds)
tFall = 0.0e-9


#Number of data points
nPoints = 50

extraPts = 10


time = np .linspace(0.0, tStart+tRise+tPulse+tFall+5e-9, nPoints+extraPts)

eField = np.zeros(nPoints+extraPts)
for i in range(nPoints):

  if (time[i] < tStart):
    eField[i] = 0.0
  elif (time[i] > tStart and time[i] <= (tStart+tRise)):
    eField[i] = (field_amplitude/tRise)*time[i]
  elif(time[i] > (tStart+tRise) and time[i] <= (tStart+tRise+tPulse)):
    eField[i] = field_amplitude
  #elif(time[i] > (tStart+tRise+tPulse)):
  #  eField[i] = (-field_amplitude/tFall)*(time[i] - (tStart+tRise+tPulse)) + field_amplitude
  else:
    eField[i] = 0.0
    
    
np.savetxt('electric_field.txt', np.array([eField, time]).T)
    
plt.plot(time, eField)
plt.show()
