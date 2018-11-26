'''
Main module for hfsolver
'''

__all__ = ['solver', 'z32']

import numpy as np
import mpmath as mp #fast, low precision implementation
from functools import partial
from scipy.optimize import brentq


_rtol = 1e-8
_z0 = np.sqrt(_rtol)/2
_g32 = partial(mp.fp.polylog, 3/2)
z32 = _g32(1)

def _fun(z, eta, nu, alpha):
    ret = _g32(z) - nu/2 + eta*np.log(z)/(2*alpha)
    return ret

def solver(nu, alpha):
    eta = np.sign(z32 - nu/2)
    if eta == 0: #avoid pushing the solver to the domain boundary
        return (z32, 0.)
    else:
        if alpha*nu <= 0.5*np.log(3*_rtol): #cut where g23(z) - z <= _rtol
            return (np.exp(alpha*nu), 0)
        else:
            z = brentq(_fun, _z0, 1, args=(eta, nu, alpha))
            xt = nu/2 - eta*np.log(z)/(2*alpha)
            xc = 0 if eta >= 0 else (nu - 2*xt)
            return (xt, xc)
