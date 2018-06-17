from math import sqrt
from scipy.constants import physical_constants, hbar, pi, atomic_mass
aB = physical_constants['Bohr radius'][0]
kB = physical_constants['Boltzmann constant'][0]

from .solver import z32

#### Sodium values
mass = 23*atomic_mass
a_scatt = 52*aB
g_int = 4*pi*hbar**2*a_scatt/mass

def lambda_therm(T):
    return sqrt((2*pi*hbar**2)/(mass*kB*T))

def n_crit(T):
    return z32/lambda_therm(T)**3

def _alpha(T):
    return g_int/(lambda_therm(T)**3 * kB*T)

def _nu(mu, T):
    return mu*lambda_therm(T)**3/g_int
