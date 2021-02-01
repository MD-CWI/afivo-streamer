import numpy as np
import matplotlib.pyplot as plt
import argparse



p = argparse.ArgumentParser(\
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument('-amplitude', type= str, required = True, help="Electric field amplitude (Volts/m)")
p.add_argument("-field_file", type=str, default='field_table.txt', \
               help="Output file name")
p.add_argument("-pulseTime", type=float, nargs = 5, required=True, help="Pulse shape times in the following order: tstart trise tpulse tfall tpost")
p.add_argument("-nPulses", type = int, default = 1, help="Number of pulses you want")
p.add_argument("-FinalTime", type=float, required = True, help="Final end time")
args = p.parse_args()

#Electric field amplitude (Volts/m)
field_amplitude = float(args.amplitude[3:])
if args.amplitude[:3] == "neg":
    field_amplitude *= -1

nPulse = args.nPulses
tFinal = args.FinalTime
def pulse(field_amplitude, times):
    tStart, tRise, tPulse, tFall, tPost = times
    nPoints = 19
    time = np.linspace(tStart, tStart+tRise+tPulse+tFall, nPoints)
    eField = np.zeros(nPoints)
    for i in range(nPoints):
      if (time[i] <= tStart):
        eField[i] = 0.0
      elif (time[i] > tStart and time[i] <= (tStart+tRise)):
        eField[i] = (field_amplitude/tRise)*(time[i] - tStart)
      elif(time[i] > (tStart+tRise) and time[i] <= (tStart+tRise+tPulse)):
        eField[i] = field_amplitude
      elif(time[i] > (tStart+tRise+tPulse) and \
           time[i] <= (tStart+tRise+tPulse+tFall)):
        eField[i] = field_amplitude + \
                    (-field_amplitude/tFall)*(time[i] - (tStart+tRise+tPulse))
    return [np.append(time,tPost), np.append(eField, 0.0)]
    
time = np.zeros(1)
eField = np.zeros(1)
tStart, tRise, tPulse, tFall, tPost = args.pulseTime
for ip in range(nPulse):
    t,e = pulse(field_amplitude, [tStart+ip*tPost, tRise, tPulse, tFall, tPost+ip*tPost])
    eField = np.append(eField, e)
    time = np.append(time, t)
eField = np.append(eField, 0.0)
time = np.append(time, tFinal)
eField = np.delete(eField,0)
time = np.delete(time,0)
head = '{0:^12s}\n{1:^14s}'.format('field_vs_time', '---------------')    
np.savetxt(args.field_file, np.array([time, eField]).T, header=head,
                                          footer='---------------', comments='')
    
plt.plot(time, eField, 'bo-')
plt.show()
