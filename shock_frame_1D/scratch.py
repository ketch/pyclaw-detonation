# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:08:52 2014

@author: maltezfaria
"""

import numpy as np
import steadyState as steadyState
import matplotlib.pylab as plt

xgrid = np.linspace(0, 100, 1000)

Ea = 1
Q = 1
gamma = 1.4


rho0, U0, p0, lam0, D, k = steadyState.steadyState(Q, Ea, gamma, xgrid)

plt.figure()

plt.plot(xgrid, rho0, xgrid, lam0, xgrid, U0+D, xgrid, p0)
plt.legend(('rho', 'lam', 'u','p'))