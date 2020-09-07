#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 01-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""

import numpy as np
from scipy.special import spence

import mpmath as mp  # fast, low precision implementation
from functools import partial

_g12 = partial(mp.fp.polylog, 1 / 2)
_g32 = partial(mp.fp.polylog, 3 / 2)
_g52 = partial(mp.fp.polylog, 5 / 2)

g12 = np.vectorize(_g12, otypes=[np.float64], doc="polylog(1/2, x). Vectorized by numpy.")
g32 = np.vectorize(_g32, otypes=[np.float64], doc="polylog(3/2, x). Vectorized by numpy.")
g52 = np.vectorize(_g52, otypes=[np.float64], doc="polylog(5/2, x). Vectorized by numpy.")


def g2(x):
    """ polylog(2, x) implemented via scipy.special.spence.
        This translation is due to scipy's definition of the dilogarithm
    """
    return spence(1.0 - x)


z2 = g2(1)
z32 = _g32(1)
z52 = _g52(1)
