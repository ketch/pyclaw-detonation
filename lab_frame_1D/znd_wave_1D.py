#!/usr/bin/env python
# encoding: utf-8
"""
Reactive euler solver
===================================

Solve the one-dimensional reactive Euler equations for inviscid, compressible flow:

.. math::
    \rho_t + (\rho u)_x & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x & = 0 \\
    (\rho E)_t + (u (\rho E + p) )_x & = 0\\
    (\rho Y)_t + (u \rho Y)_x = \rho W

The fluid is an ideal gas, with pressure given
 by :math:`p=\rho (\gamma-1) (E - u^2/2 - Q Y)` where
Q is the heat release.
"""

import numpy as np
import mpi4py.MPI as mpi
from setplot import setplot

gamma = 1.2
gamma1 = gamma - 1.
q_asympt = 1.7
theta = 2.
qheat = q_asympt*gamma1*gamma
Ea = theta/(gamma1**2)

qheat = 50
Ea = 26

if mpi.COMM_WORLD.Get_rank()==0:
    print('Heat release Q = {} and activation energy = {}'.format(qheat,Ea)) 

T_ign = 2

xmax=100
mx= 2**12
xs = 50
tfinal = 100
num_output_times = 100


def step_reaction(solver, state, dt):
    k = state.problem_data['k']
    q = state.q
    q[3,:] = q[3,:] + dt*omega(state)

def dq_step_reaction(solver, state, dt):
    k = state.problem_data['k']
    q = state.q
    dq = np.zeros(q.shape)
    dq[3,:] = dt*omega(state)

    return dq

def omega(state):
    """ Reaction rate function
    """
    k = state.problem_data['k']
    q = state.q
    rho = q[0,:]
    u   = q[1,:]/rho
    press = gamma1*(q[2,:]-0.5*q[1,:]**2/q[0,:] -
                        qheat*q[3,:])
    T = press/rho
    return -k*(q[3,:])*np.exp(Ea*(1-1/T))*(T>T_ign)    

def qinit(state,domain,xs=xs):
    import steadyState as steadyState

    grid = state.grid
    ((ixlower,ixupper),) = domain.patch._da.getRanges()


    dx = grid.delta[0] #Could be problem for non-uniform grid
    xglobal = np.linspace(dx/2,xmax-dx/2,mx)

    x =grid.x.centers
    xl = xglobal[xglobal<xs]

    rhol, Ul, pl, laml, D, k = steadyState.steadyState(qheat, Ea, gamma, xl)
    ul = Ul + D
    Yl = 1- laml

    rhor = 1.*np.ones(np.shape(xglobal[xglobal>=xs]))
    ur = 0.0*np.ones(np.shape(xglobal[xglobal>=xs]))
    pr = 1*np.ones(np.shape(xglobal[xglobal>=xs]))
    Yr = 1.*np.ones(np.shape(xglobal[xglobal>=xs]))

    rhog = np.append(rhol, rhor)
    ug = np.append(ul, ur)
    vg = 0.*ug
    pg = np.append(pl, pr)
    Yg = np.append(Yl, Yr)

    rho = rhog[ixlower:ixupper]
    u = ug[ixlower:ixupper]
    v = vg[ixlower:ixupper]
    p = pg[ixlower:ixupper]
    Y = Yg[ixlower:ixupper]

    pert = 0.000 * np.sin(2*4*np.pi*x) * (x<xs)
    state.q[0,:] = rho + pert
    state.q[1,:] = rho * u +pert
    state.q[2,:] = p/gamma1 + rho*(u**2)/2 + qheat * rho * Y  + pert
    state.q[3,:] = rho * Y

    state.problem_data['xfspeed']= D
    state.problem_data['k']= k


def setup(outdir='./_output',use_petsc=False, mx=mx, tfinal=tfinal,
          num_output_times=num_output_times, solver_type='classic'):

    from clawpack import riemann
    import reactive_euler_roe_1D

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(reactive_euler_roe_1D)
        solver.dq_src= dq_step_reaction
        solver.weno_order = 5
        #solver.lim_type = 1
        #solver.time_integrator = "SSP33"
    else:
        solver = pyclaw.ClawSolver1D(reactive_euler_roe_1D)
        solver.limiters = [1,1,1,1]
        solver.step_source = step_reaction
        solver.order = 2

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap

    solver.num_waves = 4
    solver.num_eqn = 4


    # Initialize domain
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 4
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat
    state.problem_data['Ea'] = Ea
    state.problem_data['T_ign'] = T_ign

    qinit(state,domain)

    claw = pyclaw.Controller()

    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.keep_copy = False
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    claw.setplot = setplot
    #claw.output_format = ['hdf5', 'ascii']
    claw.output_format = 'ascii'

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
