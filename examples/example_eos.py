#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 04-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""

import numpy as np
import matplotlib.pyplot as plt
from hfsolver.physics import T_crit, g_int, physics_solver_mu

n = 1e20
t = np.linspace(0, 1.6, 200)
Tc = T_crit(n)
T = t*Tc

mu, n0 = physics_solver_mu(n, T)

u = mu / g_int / n

plt.plot(t, u)
print(u)
plt.show()
