"""
potentials.py
Library of 1-D potentials V(x). Each function takes x (a number or a numpy
array) plus its own physical parameters, and returns V(x).

To study a new system: add a function here, then in config.py set
POTENTIAL to it and POTENTIAL_PARAMS to its parameters.
"""

import numpy as np


def harmonic(x, k=1.0):
    """V = 1/2 k x^2 (spring constant k; angular frequency omega = sqrt(k / m))."""
    return 0.5 * k * x**2


def double_well(x, a, b, c, d):
    """V = a x^4 + b x^3 + c x^2 + d x. a > 0 confines the particle, c < 0 makes
    two wells, b != 0 (or d != 0) makes the two wells unequal."""
    return a * x**4 + b * x**3 + c * x**2 + d * x


def flat(x):
    """V = 0. With WALLS = True this is the particle in a box."""
    return np.zeros_like(np.asarray(x, dtype=float))


def sech2(u):
    """sech^2(u) = 4 e^(-2|u|) / (1 + e^(-2|u|))^2, written so it never overflows."""
    e = np.exp(-2.0 * np.abs(u))
    return 4.0 * e / (1.0 + e)**2


def poschl_teller(x, depth, width):
    """V = -depth sech^2(x / width). Reflectionless when 2 m depth width^2 / hbar^2 = N (N + 1)
    with N integer; it then has N bound states. For hbar = m = width = 1: depth = 1 gives one
    bound state (E = -1/2), depth = 3 gives two (E = -2, -1/2). V -> 0 far away: use WALLS = True."""
    return -depth * sech2(x / width)


def dimple_in_trap(x, D, s, L):
    """Narrow well -D sech^2(x/s) at the centre of a wide, soft trap (x/L)^6."""
    return (x / L)**6 - D * sech2(x / s)


def morse(x, D, a, x0=0.0):
    """V = D (1 - exp(-a (x - x0)))^2. V -> D for x -> +infinity, so the particle
    is only bound below D: use WALLS = True or keep k_B T well below D."""
    return D * (1.0 - np.exp(-a * (x - x0)))**2
