#!/usr/bin/env python3

import argparse
import numpy as np
import sys


p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument('log_a', type=str, help='Log file')
p.add_argument('log_b', type=str, help='Log file')
p.add_argument('-atol', type=float, default=1e-8, help='Absolute tolerance')
p.add_argument('-rtol', type=float, default=1e-5, help='Relative tolerance')
p.add_argument('-skip_header', type=int, default=1,
               help='Header lines to skip')
args = p.parse_args()

log_a = np.genfromtxt(args.log_a, skip_header=args.skip_header)
log_b = np.genfromtxt(args.log_b, skip_header=args.skip_header)

if np.any(log_a.shape != log_b.shape):
    print("Log files have different shape")
    sys.exit(1)

if not np.all(np.isclose(log_a, log_b, args.rtol, args.atol)):
    print("Values in log files differ more than tolerance")
    sys.exit(1)
