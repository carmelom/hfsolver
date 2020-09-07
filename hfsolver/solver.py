'''
Main module for hfsolver
'''


import numpy as np
from scipy.optimize import brentq
from .functions import _g32, z32

__all__ = ['_solver', '_solver_mu', 'z32']

_rtol = 1e-8
_xtol = np.sqrt(3 * _rtol)

_z0 = np.sqrt(_rtol)/2


def _fun(z, eta, nu, alpha):
    ret = _g32(z) - nu / 2 + eta * np.log(z) / (2 * alpha)
    return ret


def _solver(nu, alpha):
    eta = np.sign(z32 - nu / 2)
    if eta == 0:  # avoid pushing the solver to the domain boundary
        return (z32, 0.)
    else:
        # cut where g23(z) - z <= _rtol
        if alpha * nu <= 0.5 * np.log(3 * _rtol):
            return (np.exp(alpha * nu), 0)
        else:
            z = brentq(_fun, _z0, 1, args=(eta, nu, alpha))
            xt = nu / 2 - eta * np.log(z) / (2 * alpha)
            xc = 0 if eta >= 0 else (nu - 2 * xt)
            return (xt, xc)


def _solver_mu(x, alpha):
    """
    return: (nu, xc)
    """
    # if alpha == 0:
    #     return (x, x)
    if x <= _xtol:
        return (np.log(x)/alpha, 0)
    elif x == z32:
        return (x, 0)
    else:
        if x < z32:
            F = lambda z, x: _g32(z) - x
            z = brentq(F, 0, 1, args=(x,))
            xc = 0
        else:
            F = lambda z, x, alpha: _g32(z) - x - np.log(z)/alpha
            z = brentq(F, 0, 1 - _z0, args=(x, alpha))
            xc = np.log(z)/alpha
        nu = 2*x + np.log(z)/alpha
        return nu, xc
