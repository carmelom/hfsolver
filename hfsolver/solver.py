'''
Main module for hfsolver
'''

__all__ = ['solver', 'z32']

import numpy as np
import mpmath as mp #fast, low precision implementation
from functools import partial
from scipy.optimize import brentq


_g32 = partial(mp.fp.polylog, 3/2)
z32 = _g32(1)

def _fun(z, eta, nu, alpha):
    ret =  _g32(z) - nu/2 + eta*np.log(z)/(2*alpha)
    return ret

def solver(nu, alpha):
    eta = np.sign(z32 - nu/2)
    if eta == 0:
        return (z32, 0.)
    else:
        if nu <= -1e3: #cut where numerical errors start to grow
            return (np.exp(alpha*nu), 0)
        else:
            z = brentq(_fun, 0, 1, args=(eta, nu, alpha))
            xt = nu/2 - eta*np.log(z)/(2*alpha)
            xc = 0 if eta >= 0 else (nu - 2*xt)
            return (xt, xc)
