#!/usr/bin/env python3

# Contributors:
# Jannis Teunissen
# Behnaz Bagheri

# from sympy import *
from scipy.integrate import quad, simpson
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
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
    p.add_argument('-gases', type=str, nargs='+', default=['O2'],
                   help='List of absorbing gases present')
    p.add_argument('-pressures', type=float, nargs='+', default=[0.2],
                   help='Partial pressures of gases (bar)')
    p.add_argument('-fit_range', nargs=2, type=float, default=[1e-4, 3e-3],
                   help='Distance range for fit of coefficients')
    p.add_argument('-n_modes', type=int, default=3,
                   help='Number of Helmholtz modes')
    p.add_argument('-H2O_model', type=str, choices=['Naidis', 'Aints'],
                   default='Naidis',
                   help='H2O absorption coefficients for '
                   'numerical absorption function')
    p.add_argument('-guess_amplitudes', type=float,
                   help='Initial guess for amplitude of Helmholtz modes')
    p.add_argument('-guess_lambdas', type=float,
                   help='Initial guess for lambdas of Helmholtz modes')
    p.add_argument('-fit_what', type=str, default='numerical',
                   choices=['numerical', 'Zheleznyak-H2O', 'Aints'],
                   help='What type of data/function to fit')
    p.add_argument('-fit_type', type=str, default='least_squares',
                   choices=['least_squares', 'relative', 'log'],
                   help='What type of errors to use in fit')
    p.add_argument('-ptot_for_quenching', type=float,
                   help='Total gas pressure (bar) to show quenching info')
    p.add_argument('-show_Zheleznyak', action='store_true',
                   help='Show Zheleznyak curve for air')
    p.add_argument('-show_Naidis_moist', action='store_true',
                   help='Show Naidis curve for moist air')
    p.add_argument('-show_Aints_moist', action='store_true',
                   help='Show Aints curve for moist air')
    p.add_argument('-hide_fit', action='store_true',
                   help='Hide fit in plot')
    p.add_argument('-n_points', type=int, default=300,
                   help='Number of points to use for numerical approximation')
    p.add_argument('-show_curve', metavar='A_i L_i', type=float, nargs='+',
                   help='Show curve for (A1, L1, A2, L2, ...), where A are'
                   ' amplitudes and L are lambdas (multiplied by user)')
    p.add_argument('-figure_name', type=str,
                   default='plot_of_absorption_function.png',
                   help='File name of figure')
    return p.parse_args()


args = get_args()

# Constants
ONE_TORR = 133.322368 * 1.0e-5  # bar

# Wavelength range relevant for N2-O2 mixtures
LAMBDA_MIN = 98e-9              # minimum wavelength (nm)
LAMBDA_MAX = 102.5e-9           # maximum wavelength (nm)

# Zheleznyak coefficients for air
MU_MAX = 2.0e2 / ONE_TORR       # convert to 1/(m * bar)
MU_MIN = 0.035e2 / ONE_TORR     # convert to 1/(m * bar)
pq_air_Zheleznyak = 30 * ONE_TORR  # Quenching parameter

# convert to 1/(m * bar) Naidis 2006 PSST 15 253
K_H2O = 0.26e2 / ONE_TORR
# convert to 1/(m * bar) Aints Plasma Processes and Polymers 2008 5, 672-680
K_H2O_MIN = 0.13e2 / ONE_TORR
# convert to 1/(m * bar) Aints Plasma Processes and Polymers 2008 5, 672-680
K_H2O_MAX = 0.57e2 / ONE_TORR
pq_H2O_Aints = 0.3 * ONE_TORR   # Quenching parameter


def zheleznyak_f(r):
    """Analytic absorption function at distance r (m) for dry air, unit: 1/m
    From Zheleznyak et al. 1982
    """
    return (np.exp(-MU_MIN * p['O2'] * r) -
            np.exp(-MU_MAX * p['O2'] * r)) / (r * np.log(MU_MAX / MU_MIN))


def naidis_moist_zheleznyak_f(r):
    """Analytic absorption function at distance r (m) for moist air, unit: 1/m
    This approximation is based on Naidis 2006 DOI:10.1088/0963-0252/15/2/010
    """
    nom = (np.exp(-(MU_MIN * p['O2'] + K_H2O * p['H2O']) * r) -
           np.exp(-(MU_MAX * p['O2'] + K_H2O * p['H2O']) * r))
    denom = (r * np.log(MU_MAX / MU_MIN))
    return nom/denom


def aints_moist_zheleznyak_f(r):
    """Analytic absorption function at distance r (m) for most air, unit: 1/m
    This approximation is based on Aints 2008 DOI:10.1002/ppap.200800031
    """
    nom = (np.exp(-(MU_MIN * p['O2'] + K_H2O_MIN * p['H2O']) * r) -
           np.exp(-(MU_MAX * p['O2'] + K_H2O_MAX * p['H2O']) * r))
    denom = (r * np.log((MU_MAX * p['O2'] + K_H2O_MAX * p['H2O']) /
                        (MU_MIN * p['O2'] + K_H2O_MIN * p['H2O'])))
    return nom/denom


class Absorption_O2:
    gas = "O2"
    ionizing = True

    def set_pressure(self, p):
        """Set partial pressure in bar"""
        self.p = p

    def mu(self, x):
        """Absorption function at wavelength x (nm). Unit: 1 / (m * bar)"""

        # From Zheleznyak et al. 1982
        t = (1/x - 1/LAMBDA_MAX)/(1/LAMBDA_MIN - 1/LAMBDA_MAX)
        return self.p * MU_MIN * (MU_MAX/MU_MIN)**t


class Absorption_CO2:
    """Approximation for absorption coefficient in CO2"""
    gas = "CO2"
    ionizing = False

    def set_pressure(self, p):
        """Set partial pressure in bar"""
        self.p = p

    def mu(self, x):
        """Absorption function at wavelength x (nm). Unit: 1 / (m * bar)"""
        return self.p * 1.0e2 / ONE_TORR


class Absorption_H2O:
    """Approximation for absorption coefficient in H2O"""
    gas = "H2O"
    ionizing = False

    def __init__(self, model='Naidis'):
        self.model = model
        print(self.model)

    def set_pressure(self, p):
        """Set partial pressure in bar"""
        self.p = p

    def mu(self, x):
        """Absorption function at wavelength x (nm). Unit: 1 / (m * bar)"""

        if self.model == 'Naidis':
            # This approximation is based on Naidis 2006
            # DOI:10.1088/0963-0252/15/2/010
            mu_x = self.p * 0.26e2 / ONE_TORR
        elif self.model == 'Aints':
            # This approximation is based on Aints 2008
            # DOI:10.1002/ppap.200800031
            t = (1/x - 1/LAMBDA_MAX)/(1/LAMBDA_MIN - 1/LAMBDA_MAX)
            mu_x = self.p * K_H2O_MIN * (K_H2O_MAX/K_H2O_MIN)**t
        else:
            raise ValueError('Unknown H2O absorption model')

        return mu_x


# List of all gases
all_gases = [Absorption_O2(),
             Absorption_CO2(),
             Absorption_H2O(args.H2O_model)]
all_names = [x.gas for x in all_gases]

# Construct lists of absorbing and ionizing gases present
f_absorption = []
f_ionizing = []
# Pressures
p = {'O2': 0., 'CO2': 0., 'H2O': 0.}

for gas, pressure in zip(args.gases, args.pressures):
    ix = all_names.index(gas)

    # Set partial pressure
    all_gases[ix].set_pressure(pressure)
    p[gas] = pressure

    # Append to list of absorption functions
    f_absorption.append(all_gases[ix].mu)

    # If ionizing, add to list of ionizing gases
    if all_gases[ix].ionizing:
        f_ionizing.append(all_gases[ix].mu)


def integrand(x, r):
    """Compute integrand that can be integrated over wavelength

    :param x: wavelength (nm)
    :param r: distance (m)

    """
    ionization = sum([f(x) for f in f_ionizing])
    absorption = sum([f(x) for f in f_absorption])
    return ionization * np.exp(-absorption * r)


r = np.linspace(args.fit_range[0], args.fit_range[1], args.n_points)
f_numerical = np.zeros(args.n_points)

# Construct absorption function numerically integrated over wavelength
for i in range(args.n_points):
    tmp = quad(integrand, LAMBDA_MIN, LAMBDA_MAX, args=(r[i]))
    f_numerical[i] = tmp[0] / (LAMBDA_MAX - LAMBDA_MIN)


if args.fit_what == 'numerical':
    f_to_fit = f_numerical
elif args.fit_what == 'Zheleznyak-H2O':
    f_to_fit = naidis_moist_zheleznyak_f(r)
elif args.fit_what == 'Aints':
    f_to_fit = aints_moist_zheleznyak_f(r)
    # Show quenching information
    if args.ptot_for_quenching is not None:
        Q = (1 + (args.ptot_for_quenching - p['H2O'])/pq_air_Zheleznyak +
             p['H2O']/pq_H2O_Aints)**-1
        pq_effective = -Q * args.ptot_for_quenching/(Q - 1)
        print('Effective quenching pressure for Aints model: '
              f'{pq_effective:.5e} bar')
else:
    raise ValueError('Invalid argument for fit_what')


lambda_guess = -np.log(f_to_fit[-2]/f_to_fit[-1])/(r[-2] - r[-1])
amplitude_guess = lambda_guess**2 / args.n_modes

# Construct interpolating function
f_func = interp1d(r, f_to_fit)

# Guesses for the fit
fit_guess = np.ones((2 * args.n_modes))
if args.guess_amplitudes is not None:
    fit_guess[0::2] = args.guess_amplitudes
else:
    fit_guess[0::2] = amplitude_guess

if args.guess_lambdas is not None:
    fit_guess[1::2] = args.guess_lambdas
else:
    fit_guess[1::2] = lambda_guess


# Function for fitting with a variable number of exponential functions
def fit_func(x, *fitargs):
    amplitudes = fitargs[0::2]
    lambdas = fitargs[1::2]
    val = np.zeros_like(x)
    for c1, c2 in zip(amplitudes, lambdas):
        val += x * c1 * np.exp(-c2 * x)
    return val


def fit_func_relative(x, *fitargs):
    return fit_func(x, *fitargs)/f_func(x)


def fit_func_log(x, *fitargs):
    return np.log(fit_func(x, *fitargs))


try:
    if args.fit_type == 'relative':
        popt, pcov = curve_fit(fit_func_relative, r,
                               np.ones(args.n_points), p0=fit_guess)
    if args.fit_type == 'log':
        popt, pcov = curve_fit(fit_func_log, r,
                               np.log(f_to_fit), p0=fit_guess)
    else:
        popt, pcov = curve_fit(fit_func, r, f_to_fit, p0=fit_guess)
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

print('Fit sigma:     {:.5e}'.format(
    np.mean(np.sqrt(np.diag(pcov)) / np.abs(popt))))
print('Fit range (m): {:.5e} -- {:.5e}'.format(
    args.fit_range[0], args.fit_range[1]))

f_to_display = []
if args.show_Zheleznyak:
    f_to_display.append(['Zheleznyak air', zheleznyak_f])
if args.show_Naidis_moist:
    f_to_display.append(['Naidis moist air', naidis_moist_zheleznyak_f])
if args.show_Aints_moist:
    f_to_display.append(['Aints moist air', aints_moist_zheleznyak_f])

print('Integrals of absorption functions over fit range:')
print(f'{"Numerical":<20} {simpson(f_numerical, r):12.5f}')
for name, f in f_to_display:
    print(f'{name:<20} {simpson(f(r), r):12.5f}')

fig, ax = plt.subplots(1, 2, layout='constrained', figsize=(7, 4))

ax[0].set_xlabel('r (m)')
ax[0].set_ylabel('absorption function (1/m)')
ax[0].set_title('Logarithmic scale')
ax[0].semilogy(r, f_numerical, '.-', label='numerical')

if not args.hide_fit:
    ax[0].semilogy(r, fit_func(r, *popt),
                   label='fit ({}-term)'.format(args.n_modes))

for name, f in f_to_display:
    ax[0].semilogy(r, f(r), '--', label=name)

if args.show_curve:
    ax[0].semilogy(r, fit_func(r, *args.show_curve), label='custom')
ax[0].legend()

ax[1].set_xlabel('r (m)')
ax[1].set_ylabel('absorption function (1/m)')
ax[1].set_title('Linear scale')
ax[1].plot(r, f_numerical, '.-', label='numerical')

if not args.hide_fit:
    ax[1].plot(r, fit_func(r, *popt),
               label='fit ({}-term)'.format(args.n_modes))

for name, f in f_to_display:
    ax[1].plot(r, f(r), '--', label=name)

if args.show_curve:
    ax[1].plot(r, fit_func(r, *args.show_curve), label='custom')
ax[1].legend()

plt.savefig(args.figure_name, bbox_inches='tight', dpi=200)
print('Saved {}'.format(args.figure_name))
