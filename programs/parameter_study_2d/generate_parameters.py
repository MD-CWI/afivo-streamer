#!/usr/bin/env python3

import numpy as np
from numpy.random import default_rng

rng = default_rng()

# Start of user code
n_samples = 1000
cfg_file = "malagon_dataset.cfg"

parameters = {
    'field_amplitude': rng.uniform(1.5e6, 2.5e6, n_samples),
    'field_rod_r1': np.vstack(
        [np.zeros(n_samples),
         rng.uniform(0.125, 0.175, n_samples)]).T
    }

# Generate commands
names = parameters.keys()
commands = []

for index in range(n_samples):
    values = [parameters[name][index] for name in names]

    # Convert array arguments
    for i, v in enumerate(values):
        if isinstance(v, (list, tuple, np.ndarray)):
            values[i] = ' '.join(map(str, v))

    arguments = [f"-{a}='{b}'" for a, b in zip(names, values)]
    arguments.append(f'-output%name+=_run{index}')
    commands.append(f'./streamer {cfg_file} ' + ' '.join(arguments))

with open('commands.txt', 'w') as f:
    for c in commands:
        f.write(c + '\n')
