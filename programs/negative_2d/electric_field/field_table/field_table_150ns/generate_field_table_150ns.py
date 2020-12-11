import numpy as np
import matplotlib.pyplot as plt
import argparse


p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument('-field_file', type=str, default='field_table.txt', help='Output file name')
args = p.parse_args()


# Initial field before any decay (V/m)
initial_field=2.5e6

# Minimal field (V/m)
min_field=0.8e6

# Decay time (s)
decay_time = 10e-9

# Decay slope
decay_slope = -0.9e14

# Pulse duration (s)
tPulse = 23.5e-9

# Fall time (s)
tFall = 126.5e-9

# Total time (s)
tTotal = tPulse + tFall

# Number of data points
nPoints = 2000

time = np.linspace(0.0, tTotal, nPoints)

t=np.maximum(time-tPulse, 0.0)


# Electric field profile
decay_field_profile={
  'constant': 0,  # Constant electric field versus time  
  'exponential': 1,  # Exponential decay electric field  versus time
  'linear_across_zero': 2,  # Linear decay electric field versus time (across zero)  
  'linear_over_zero': 3, # Linear decay electric field versus time (over zero)
  'step': 4  # Step decay electric field versus time
 } 
profile=decay_field_profile.get('exponential')

if (profile==1):
  eField=min_field+(initial_field-min_field)*np.exp(-t/decay_time)
elif (profile==2):
  eField=initial_field+decay_slope*t
  for i in range(nPoints):
    if (eField[i]<=-initial_field):
      eField[i] = -initial_field
elif (profile==3):
  eField=initial_field+decay_slope*t
  for i in range(nPoints):
    if (eField[i]<=0.0):
      eField[i] = 0.0
elif (profile==4):
  eField = np.zeros(nPoints)
  for i in range(nPoints):
    if (time[i]<=tPulse):
      eField[i] = initial_field
    else:
      eField[i]=min_field
else:
  eField = np.zeros(nPoints)
  for i in range(nPoints):
    eField[i]=initial_field


head = '{0:^12s}\n{1:^14s}'.format('field_vs_time', '---------------')    
np.savetxt(args.field_file, np.array([time, eField]).T, header=head, footer='---------------', comments='')


# Convert to kV/cm
eField=eField/1e5

# Convert to ns
time=time/1e-9

plt.figure(figsize=(10,6))
plt.plot(time, eField, label='Electric field', color='red', linestyle='-', linewidth=2, marker='D', markersize=0)
plt.tick_params(labelsize=14)
plt.xlabel('Time (ns)', fontsize=14)
plt.ylabel('Electric field (kV/cm)', fontsize=14)
plt.grid(True)
plt.grid(linestyle='-.')
plt.legend(bbox_to_anchor=(0.95,0.95), frameon=False, fontsize=16 )
plt.savefig('field_table.jpg')
plt.show()


