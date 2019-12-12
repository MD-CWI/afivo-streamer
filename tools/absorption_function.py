#!/usr/bin/env python3

# Author: Jannis Teunissen, CWI

# from sympy import *
from scipy.integrate import quad, simps
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys


def get_args():
    # Get and parse the command line arguments
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Compute absorption coefficients for
        Helmholtz photoionization''')
    p.add_argument('-p_O2', type=float, default=0.2,
                   help='Partial pressure of O2 (bar)')
    p.add_argument('-p_CO2', type=float, default=0.0,
                   help='Partial pressure of CO2 (bar)')
    p.add_argument('-fit_range', nargs=2, type=float, default=[1e-4, 3e-3],
                   help='Distance range for fit of coefficients')
    p.add_argument('-n_modes', type=int, default=2,
                   help='Number of Helmholtz modes')
    p.add_argument('-guess_amplitudes', type=float, default=1e6,
                   help='Initial guess for amplitude of Helmholtz modes')
    p.add_argument('-guess_lambdas', type=float, default=1e3,
                   help='Initial guess for lambdas of Helmholtz modes')
    p.add_argument('-n_points', type=int, default=300,
                   help='Number of points to use for numerical approximation')
    return p.parse_args()


args = get_args()
one_torr = 133.322368 * 1.0e-5  # bar
lambda_min = 98e-9              # nm
lambda_max = 102.5e-9           # nm
mu_max = 2.0e2 / one_torr       # convert to 1/(m * bar)
mu_min = 0.035e2 / one_torr     # convert to 1/(m * bar)


# Zheleznyak approximation for absorption coefficient
# Unit: 1 / (m * bar)
def mu_O2(x):
    t = (1/x - 1/lambda_max)/(1/lambda_min - 1/lambda_max)
    return mu_min * (mu_max/mu_min)**t


# Approximation for absorption coefficient in CO2
# Unit: 1 / (m * bar)
def mu_CO2(x):
    return 1.0e2 / one_torr


def mu_all(x):
    return mu_O2(x) * args.p_O2 + mu_CO2(x) * args.p_CO2


def nominator(x, r):
    return mu_O2(x) * args.p_O2 * np.exp(-mu_all(x) * r)


def zheleznyak_f(r):
    return (np.exp(-mu_min * args.p_O2 * r) -
            np.exp(-mu_max * args.p_O2 * r)) / (r * np.log(mu_max / mu_min))


r = np.linspace(args.fit_range[0], args.fit_range[1], args.n_points)
g = zheleznyak_f(r)
f = np.zeros(args.n_points)

for i in range(args.n_points):
    f[i] = quad(nominator, lambda_min, lambda_max,
                args=(r[i]))[0] / \
        (lambda_max - lambda_min)


# Function for fitting with a variable number of exponential functions
def fit_func(x, *args):
    amplitudes = args[0::2]
    lambdas = args[1::2]
    val = np.zeros_like(x)
    for c1, c2 in zip(amplitudes, lambdas):
        val += x * c1 * np.exp(-c2 * x)
    return val


fit_guess = np.ones((2 * args.n_modes))
fit_guess[0::2] = args.guess_amplitudes
fit_guess[1::2] = args.guess_lambdas

try:
    popt, pcov = curve_fit(fit_func, r, f, p0=fit_guess)
except RuntimeError as e:
    print('No convergence, adjust guess_amplitudes and/or guess_lambdas')
    print(e)
    sys.exit(1)

print('{:>15s} {:>15s}'.format('amplitude', 'lambda'))
print('----------------------------------------')
amplitudes = popt[0::2]
lambdas = popt[1::2]

for c1, c2 in zip(amplitudes, lambdas):
    print('{:15.5e} {:15.5e}'.format(c1, c2))
print('----------------------------------------')

# print('Fit sigma {:.5e}'.format(
#     np.mean(np.sqrt(np.diag(pcov)) / np.abs(popt))))
print('Fit range (m): {:.5e} -- {:.5e}'.format(
    args.fit_range[0], args.fit_range[1]))
print('Integral of absorption function over fit range: {:.5e}'.format(
    simps(f, r)))

fig, ax = plt.subplots(1, 2)

plt.subplot(121)
plt.xlabel('r (m)')
plt.ylabel('absorption function (1/m)')
plt.title('Logarithmic scale')
plt.semilogy(r, f, '.', label='numerical')
plt.semilogy(r, g, label='Zheleznyak air')
plt.semilogy(r, fit_func(r, *popt), label='fit ({}-term)'.format(args.n_modes))
plt.legend()
plt.subplot(122)
plt.xlabel('r (m)')
plt.ylabel('absorption function (1/m)')
plt.title('Linear scale')
plt.plot(r, f, '.', label='numerical')
plt.plot(r, g, label='Zheleznyak air')
plt.plot(r, fit_func(r, *popt), label='fit ({}-term)'.format(args.n_modes))
plt.legend()

fig.tight_layout()
fname = 'plot_of_absorption_function.png'
plt.savefig(fname, bbox_inches='tight', dpi=200)
print('Saved {}'.format(fname))
