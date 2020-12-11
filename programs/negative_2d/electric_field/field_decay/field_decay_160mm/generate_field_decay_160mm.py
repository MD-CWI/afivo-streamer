import numpy as np
import matplotlib.pyplot as plt
import argparse


p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument('-field_file', type=str, default='2.1e6_0.6e6_50e-3_145e-3.txt', help='Output file name')
args = p.parse_args()


# Initial field before any decay (V/m)
initial_field=2.1e6

# Minimal field (V/m)
min_field=0.6e6

# Decay distance (m)
decay_distance = 50e-3

# Decay slope
decay_slope = -0.5e8

# Pulse length (m)
zPulse = 15e-3

# Fall length (m)
zFall = 145e-3

# Total length (m)
zTotal = zPulse + zFall

# Number of data points
nPoints =800

length = np.linspace(0.0, zTotal, nPoints)

L=np.maximum(length-zPulse, 0.0)


# Electric field profile
decay_field_profile={
  'constant': 0,  # Constant electric field versus length  
  'exponential': 1,  # Exponential decay electric field  versus length
  'linear_across_zero': 2,  # Linear decay electric field versus length (across zero)  
  'linear_over_zero': 3, # Linear decay electric field versus length (over zero)
  'step': 4  # Step decay electric field versus length
 } 
profile=decay_field_profile.get('exponential')

if (profile==1):
  eField=min_field+(initial_field-min_field)*np.exp(-L/decay_distance)
elif (profile==2):
  eField=initial_field+decay_slope*L
  for i in range(nPoints):
    if (eField[i]<=-initial_field):
      eField[i] = -initial_field
elif (profile==3):
  eField=initial_field+decay_slope*L
  for i in range(nPoints):
    if (eField[i]<=0.0):
      eField[i] = 0.0
elif (profile==4):
  eField = np.zeros(nPoints)
  for i in range(nPoints):
    if (length[i]<=zPulse):
      eField[i] = initial_field
    else:
      eField[i]=min_field
else:
  eField = np.zeros(nPoints)
  for i in range(nPoints):
    eField[i]=initial_field


head = '{0:^12s}\n{1:^14s}'.format('field_vs_length', '---------------')    
np.savetxt(args.field_file, np.array([0.16-length, eField]).T, header=head, footer='---------------', comments='')


# Convert to kV/cm
eField=eField/1e5

# Convert to mm
length=length/1e-3

plt.figure(figsize=(10,6))
plt.plot(160-length, eField, label='2.1e6_0.6e6_50e-3_145e-3', color='red', linestyle='-', linewidth=2, marker='D', markersize=0)
plt.tick_params(labelsize=14)
plt.xlim((170,-10))
plt.xlabel('z-axis (mm)', fontsize=14)
plt.ylabel('Electric field (kV/cm)', fontsize=14)
plt.grid(True)
plt.grid(linestyle='-.')
plt.legend(bbox_to_anchor=(0.95,0.95), frameon=False, fontsize=16 )
plt.savefig('2.1e6_0.6e6_50e-3_145e-3.jpg')
plt.show()


