#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Script to generate commands for a sensitivity study
    Authors: Hemaditya Malla, Jannis Teunissen''')
parser.add_argument('cfg_file', type=str, help='Base config file')
parser.add_argument('-command_file', type=str, default='commands.txt',
                    help='Base config file')
parser.add_argument('-ix_range', type=int, nargs=2, required=True,
                    help='Index range of reactions to modify')
parser.add_argument('-rate_factors', type=float, nargs='+',
                    default=[0.8, 1.2],
                    help='List of reaction rate factors')
args = parser.parse_args()

commands = []

# Add base case (with modified output name)
arguments = [f'-output%name+=_ix{0:04d}_fac{1.0}']
commands.append(f'./streamer {args.cfg_file} ' + ' '.join(arguments))

for index in range(args.ix_range[0], args.ix_range[1]+1):
    for fac in args.rate_factors:
        arguments = [f'-input_data%modified_reaction_ix={index}',
                     f'-input_data%modified_rate_factors={fac}',
                     f'-output%name+=_ix{index:04d}_fac{fac}']
        commands.append(f'./streamer {args.cfg_file} ' + ' '.join(arguments))

with open(args.command_file, 'w') as f:
    for c in commands:
        f.write(c + '\n')
