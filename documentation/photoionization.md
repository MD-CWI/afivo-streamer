# Photoionization

[TOC]

# Introduction

Photoionization can play an important role in electrical discharges in air and
other gases. Two approaches for photoionization are included in afivo-streamer: the Helmholtz approximation and a Monte Carlo procedure to generate and absorb photons. These approaches are briefly explained below, and in more detail in the following references:

* \cite Luque_2007 \cite Bourdon_2007 Introduction of the Helmholtz approach
* \cite Bagheri_2018a Description of Helmholtz coefficients for air
* \cite Chanrion_2008 \cite Teunissen_2016 - Description of the Monte Carlo approximation to the Zheleznyak model (and tho photoionization in general)
* \cite Zheleznyak_1982 Description of Zheleznyak approximation for photoionization in air
* \cite Bagheri_2019 Comparison between Helmholtz and Monte Carlo photoionization

# General parameters

There are several parameters that control the production of photoionization, which are valid for both the Helmholtz and Monte Carlo approach:

\snippet m_photoi.f90 photoi_general_parameters

# Helmholtz approach

The basic idea is to approximate the absorption function for photons with a series expansion. The user can either use one of the built-in expansions, or define custom coefficients through the following parameters:

\snippet m_photoi_helmh.f90 helmholtz_parameters

For further details, see the source code of @ref m_photoi_helmh::photoi_helmh_initialize

There is a tool available to generate Helmholtz coefficients for a new gas, see @ref documentation/tools.md

# Monte Carlo approach

The basic idea of the Monte Carlo approach is to use random numbers to sample the production of photons as well as their absorption. The main parameters that control Monte Carlo photoionization are:

\snippet m_photoi_mc.f90 photoi_mc_parameters

These parameters are somewhat difficult to understand, but for almost all applications the following (default) parameters should work:

    [photoi_mc]
        # At which grid spacing photons are absorbed compared to their mean distance:
        absorp_fac = 1e-9

        # Whether a constant grid spacing is used for photoionization:
        const_dx = T

        # Minimum grid spacing for photoionization:
        min_dx =  1.0000E-09

        # Minimal photon weight (default: 1.0):
        min_weight =  1.0000E+00

        # Maximum number of discrete photons to use:
        num_photons = 5000000

        # Whether physical photons are used:
        physical_photons = T

This will cause photons to be absorbed on the finest grid, and in case of axisymmetric or 3D simulations, the produced photons will have a weight of 1.0 (corresponding to physical photons).
