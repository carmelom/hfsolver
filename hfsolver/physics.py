import numpy as np
from scipy.constants import physical_constants, hbar, pi, atomic_mass
aB = physical_constants['Bohr radius'][0]
kB = physical_constants['Boltzmann constant'][0]

from .solver import solver, z32
from scipy.interpolate import interp1d
from scipy.integrate import trapz

#### Sodium values
mass = 23*atomic_mass
a_scatt = 52*aB
g_int = 4*pi*hbar**2*a_scatt/mass

def lambda_therm(T):
    return np.sqrt((2*pi*hbar**2)/(mass*kB*T))

def n_crit(T):
    return z32/lambda_therm(T)**3

def T_crit(n):
    return 2*pi*hbar**2 * (n/z32)**(2./3)/(mass*kB)

def _alpha(T):
    return g_int/(lambda_therm(T)**3 * kB*T)

def _nu(mu, T):
    return mu*lambda_therm(T)**3/g_int

def radial_harmonic_trap(r, omega_ho=2*pi):
    return 0.5 * mass * omega_ho**2 * r**2

def integrate_N_harmonic_trap(mu0, omega_ho, n_radial, mu_radial):
    return -4*pi/np.sqrt(mass*omega_ho**2)**3 * trapz(n_radial * np.sqrt(2*(mu0 - mu_radial)), mu_radial)



def solver_harmonic_trap(mu0, T, omega_ho, r=None, dr=1e-6, Rmax=None):
    """Returns interpolating functions to calculate the density as a function of ``mu(r)``.
    Those can be evaluated on an arbitrary grid of values for ``mu``.

    All physical values are expected to be given in SI units.

    Parameters
    ----------
    mu0 : float
        Chemical potential at trap center.
    T : float
        Temperature.
    omega_ho : float
        Averaged trap frequency.
    dr : float, optional
        Spatial resolution on which to compute the interpolating grid for ``mu``.
        Defaults to 1 um.

    Returns
    -------
    fun_n : ``scipy.interpolate.interp1d`` object
        Interpolating function that returns n (the total density) as a function of mu
    fun_n0 : ``scipy.interpolate.interp1d`` object
        Interpolating function that returns n0 (BEC density) as a function of mu
    mu : ndarray
        The grid of interpolating points for mu
    r : ndarray
        The (radial, evenly spaced) spatial grid of points corresponding to ``mu``
    alpha : float
        Temperature parameter
    """
    lt = lambda_therm(T)
    alpha = _alpha(T)
    if r is None:
        if Rmax is None:
            if mu0 > 0:
                Rmax = np.sqrt(30*mu0/mass)/omega_ho # almost 4\sigma
            else:
                Rmax = 3*np.sqrt(kB*T/mass)/omega_ho # 4\sigma thermal
        else:
            Rmax = Rmax
        r = np.arange(0, Rmax, dr)
    else:
        r = r
    mu = mu0 - radial_harmonic_trap(r, omega_ho)
    nu = _nu(mu, T)
    x_t = np.empty(nu.shape)
    x_c = np.empty(nu.shape)
    for j, n in enumerate(nu):
        x_t[j], x_c[j] = solver(n, alpha)
    x_tot = x_t + x_c
    _n = x_tot.astype(float)/lt**3
    _n0 = x_c.astype(float)/lt**3
    fun_n = interp1d(mu, _n, kind='cubic')
    fun_n0 = interp1d(mu, _n0, kind='cubic')
    return fun_n, fun_n0, mu, r, alpha
