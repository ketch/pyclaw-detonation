# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate
import matplotlib.pylab as plt

def steadyState(Q, E, gamma, xgrid):
    """Computes the steady state given the heat release and activation energy
    Ea is the activation energy
    Q is the heat release
    """
    global b, k, Ea

    x = xgrid[::-1]
    q = (gamma**2-1.)*Q/2.
    D = np.sqrt(gamma+q)+np.sqrt(q)
    b = (D**2-gamma)/(gamma*(1.+D**2))
    k = 1
    Ea = E

    def omega(lam, T):
        global Ea, k

        return k*np.abs((1-lam))*np.exp(-Ea/T)*(lam<=1)

    def dxdlam(lam):
        global Ea, b, k

        delta = b*np.sqrt(1-lam)
        rho = (gamma+1)*D**2/(gamma*(1+D**2)*(1-delta))
        p = (1+D**2)*(1+gamma*delta)/(gamma+1)
        U = -gamma*(1+D**2)*(1-delta)/((gamma+1)*D)
        T = p/rho
        return U/omega(lam, T)

    k, err = integrate.quad(dxdlam, 0, 0.5)
    k = np.abs(k)

#    print("Half reaction lenght set to {}".format(k))
#    print("Frame speed set to {} ".format(D))


    def dlamdx(lam, x):
        global Ea, b, k

        delta = b*np.sqrt((1-np.minimum(lam,1)))
        rho = (gamma+1)*D**2/(gamma*(1+D**2)*(1-delta))
        p = (1+D**2)*(1+gamma*delta)/(gamma+1)
        U = -gamma*(1+D**2)*(1-delta)/((gamma+1)*D)
        T = p/rho
        return omega(lam, T)/U

    sol = integrate.odeint(dlamdx, 0, x,rtol=10**-13,atol=10**-13)
    lam = sol[::-1, 0]

    whereAreNaNs = np.isnan(lam);
    lam[whereAreNaNs] = 1;

#    lam = lam*(lam < 1) + (lam > 1)

    delta = b*np.sqrt(np.abs(1-lam))*(lam<=1)
    rho = (gamma+1)*D**2/(gamma*(1+D**2)*(1-delta))
    p = (1+D**2)*(1+gamma*delta)/(gamma+1)
    U = -gamma*(1+D**2)*(1-delta)/((gamma+1)*D)

    return rho, U, p, lam, D, k
