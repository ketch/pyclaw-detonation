# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 14:11:55 2014

@author: maltezfaria
"""

import numpy as np
import time
import matplotlib.pylab as plt
from clawpack import pyclaw
#import clawpack.petclaw as pyclaw

from clawpack import riemann
import reactive_euler_efix_roe_1D as rp_solver
import steadyState as steadyState

plt.close('all')

def step_reaction(solver, state, dt):
    q = state.q    
    q[3,:] = q[3,:] + dt*omega(state)

def omega(state):
    """ Reaction rate function
    """
    q = state.q
    k = state.problem_data['k']
    gamma = state.problem_data['gamma']
    qheat = state.problem_data['qheat']    
    Ea = state.problem_data['Ea'] 
    T_ign = state.problem_data['T_ign']         
    gamma1 = gamma-1
    pressure = gamma1*(q[2,:]-0.5*q[1,:]**2/q[0,:] -
                        qheat*q[3,:])
    T = pressure / q[0,:]
    return -k*(q[3,:])*np.exp(Ea*(1-1/T))*(T>T_ign)

def main_solver_laf_ree(Ea = 20, qheat = 50, gamma = 1.2, T_ign = 1.5, 
                        mx=1000, xmax=30, tmax=10,dtout=1, 
                        outdir='./_output', iplot=1):
                            
    gamma1 = gamma-1
    # Set up the solver object
    solver = pyclaw.ClawSolver1D(rp_solver)    
    solver.step_source = step_reaction
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    # solver.user_bc_upper = custom_bc
    solver.num_waves = 4
    solver.num_eqn = 4

    # Set the domain
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    domain = pyclaw.Domain([x])

    # Set state    
    state = pyclaw.State(domain,solver.num_eqn)

    # Set initial condition
    x =state.grid.x.centers
    xs = xmax-5
    xdet = x[x<xs]
    rhol, Ul, pl, laml, D, k = steadyState.steadyState(qheat, Ea, gamma, xdet)
    ul = Ul + D
    Yl = 1-laml
    rhor = 1.*np.ones(np.shape(x[x>=xs]))
    ur = 0.0*np.ones(np.shape(x[x>=xs]))
    pr = 1*np.ones(np.shape(x[x>=xs]))
    Yr = 1.*np.ones(np.shape(x[x>=xs]))
    rho = np.append(rhol, rhor)
    u = np.append(ul, ur)
    p = np.append(pl, pr)
    Y = np.append(Yl, Yr)

    plt.plot(x,u)
    
    state.q[0,:] = rho
    state.q[1,:] = rho * u
    state.q[2,:] = p/gamma1 + rho*u**2/2 + qheat * rho * Y
    state.q[3,:] = rho * Y
    state.mF = 1

    # Fill dictory of problem data
    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat
    state.problem_data['Ea'] = Ea
    state.problem_data['T_ign'] = T_ign
    state.problem_data['xfspeed'] = D
    state.problem_data['k'] = k
    
    # Set up controller   

    claw = pyclaw.Controller()
    claw.tfinal = tmax
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    #Set up data writting and plotting 
    claw.output_style = 1
    nout = np.ceil(tmax/dtout)
    claw.num_output_times = nout    
    claw.outdir = outdir
    claw.keep_copy = iplot    
    
    #Run the simulation
    claw.run()    
    
    #Save QOF
    np.save('qof',QOF)
    
    return QOF

    
xmax=50
plt.ion()
ax = plt.gca()    
h = ax.plot([],[],'ro-')
plt.xlim([0, xmax])
plt.ylim([0, 2])

gamma = 1.2
gamma1 = gamma-1
qheat = 50
Ea = 26

tmax = 5
mx = 1000

main_solver_laf_ree(Ea = Ea, qheat = qheat, gamma = gamma, T_ign = 1.05, 
                        mx=mx, xmax=xmax, tmax=tmax, dtout=1, 
                        outdir='./_output', iplot=1)   

