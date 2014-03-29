# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:08:52 2014

@author: maltezfaria
"""

import numpy as np
import steadyState as steadyState
import matplotlib.pylab as plt

xgrid = np.linspace(0, 30, 1000)

gamma = 1.2
Ea = 10
Q = 10*(gamma/(gamma-1))

rho0, U0, p0, lam0, D, k = steadyState.steadyState(Q, Ea, gamma, xgrid)

plt.figure()

plt.plot(xgrid, rho0, xgrid, lam0, xgrid, U0+D, xgrid, p0)
plt.legend(('rho', 'lam', 'u','p'))
plt.show()


plt.figure()
plt.plot(xgrid, p0/rho0)
plt.legend('T')
plt.show()

#print(p0[-1]/rho0[-1])

plt.figure()
plt.plot(xgrid, U0 + np.sqrt(gamma*p0/rho0))
plt.legend('u+c-D', loc=2)
plt.legend
plt.show()