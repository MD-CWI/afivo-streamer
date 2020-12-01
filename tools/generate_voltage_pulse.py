import numpy as np
import matplotlib.pyplot as plt
import argparse



p = argparse.ArgumentParser(\
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#p.add_argument('-amplitude', type= float, help="Electric field amplitude (Volts/m)")
p.add_argument("-field_file", type=str, default='field_table.txt', \
               help="Output file name")
args = p.parse_args()

#Electric field amplitude (Volts/m)
field_amplitude = -0.18e+7  #-0.1975308e+7

#Pulse start time
tStart = 115.0e-9

#Pulse duration (seconds) (amplitude duration)
tPulse = 10.0e-9


#Rise time (seconds)
tRise = 2.5e-9

#Fall time (seconds)
tFall = 2.5e-9

#Post pulse time (seconds)
tPost = 100.0e-9


#Total time
tTotal = tStart + tRise + tPulse + tFall + tPost


#Number of data points
nPoints = 100



time = np.linspace(0.0, tTotal, nPoints)

eField = np.zeros(nPoints)
for i in range(nPoints):
  if (time[i] <= tStart):
    #print("test")
    eField[i] = 0.0
  elif (time[i] > tStart and time[i] <= (tStart+tRise)):
    #print("test")
    eField[i] = (field_amplitude/tRise)*(time[i] - tStart)
  elif(time[i] > (tStart+tRise) and time[i] <= (tStart+tRise+tPulse)):
    #print("test")
    eField[i] = field_amplitude
  elif(time[i] > (tStart+tRise+tPulse) and \
       time[i] <= (tStart+tRise+tPulse+tFall)):
    #print("test")
    eField[i] = field_amplitude + \
                (-field_amplitude/tFall)*(time[i] - (tStart+tRise+tPulse))
  else:
    #print("test")
    eField[i] = 0.0
    
head = '{0:^12s}\n{1:^14s}'.format('field_vs_time', '---------------')    
np.savetxt(args.field_file, np.array([time, eField]).T, header=head,
                                          footer='---------------', comments='')
    
plt.plot(time, eField)
plt.show()
