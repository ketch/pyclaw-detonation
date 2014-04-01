#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from setplot import setplot

gamma = 1.2
gamma1 = gamma - 1.
q_asympt = 1.7
theta = 2.
qheat = q_asympt*gamma1*gamma
Ea = theta/(gamma1**2)

#qheat = 50
#Ea = 20

print('Heat release Q={} and activation energy={}'.format(qheat,Ea))

T_ign = 1.005

xmax = 40.
ymax = 100.
xs = 25.

def b4step(solver, state):
    pass

def qinit(state,xs=xs):
    import steadyState as steadyState

    grid = state.grid
    
    # Create an array with fortran native ordering
    x =grid.x.centers
    y =grid.y.centers
    xl = x[x<xs]

    rhol, Ul, pl, laml, D, k = steadyState.steadyState(qheat, Ea, gamma, xl)
    ul = Ul + D
    Yl = 1- laml

    rhor = 1.*np.ones(np.shape(x[x>=xs]))
    ur = 0.0*np.ones(np.shape(x[x>=xs]))
    pr = 1*np.ones(np.shape(x[x>=xs]))
    Yr = 1.*np.ones(np.shape(x[x>=xs]))

    rho = np.append(rhol, rhor)
    u = np.append(ul, ur)
    v = 0.*u
    p = np.append(pl, pr)
    Y = np.append(Yl, Yr)

    for i in range(y.size):
        pert = 0.00001 * np.sin(2*np.pi*y[i]/ymax)
        state.q[0,:,i] = rho + pert
        state.q[1,:,i] = rho * u + pert
        state.q[2,:,i] = 0. + pert
        state.q[3,:,i] = p/gamma1 + rho*(u**2 + v**2)/2 + qheat * rho * Y + pert
        state.q[4,:,i] = rho * Y 
        
    state.problem_data['fspeed']= D
    state.problem_data['k']= k
    
def step_Euler_reaction(solver,state,dt):
    """
    Source terms for reactive Euler equations.
    This is a Clawpack-style source term routine.
    """
    q = state.q
    q[4,:,:] = q[4,:,:] + dt*omega(state)

def omega(state):
    """ Reaction rate function
    """
    k = state.problem_data['k']
    q = state.q
    rho = q[0,:,:]
    u   = q[1,:,:]/rho
    v   = q[2,:,:]/rho
    press = gamma1 *  (q[3,:,:] - 0.5*rho*(u**2 + v**2) - qheat*q[4,:,:])
    T = press/rho
    
    return -k*(q[4,:,:])*np.exp(Ea*(1-1/T))*(T>T_ign)
    

def setup(use_petsc=False,kernel_language='Fortran',solver_type='classic',
          outdir='_output', disable_output=False, mx=200, my=1000, tfinal=500,
          num_output_times = 500):
    """
    Solve the reactive Euler equations of compressible fluid dynamics.
    """
    from clawpack import riemann

    import reactive_euler_roe_2D

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    
    solver = pyclaw.ClawSolver2D(reactive_euler_roe_2D)
    solver.limiters = [4,4,4,4,2]
    solver.step_source=step_Euler_reaction
    solver.before_step = b4step

    # Initialize domain
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    y = pyclaw.Dimension('y',0.0,ymax,my)
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,solver.num_eqn)
    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat

    # The k and fspeed problem data are set inside qinit
    qinit(state)

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.dt_initial=0.005
    solver.source_split = 1
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.periodic
    solver.bc_upper[1]=pyclaw.BC.periodic

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    #claw.output_format = 'netcdf'

    claw.keep_copy = True

    if disable_output:
        claw.output_format = None

    claw.tfinal = tfinal
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    claw.setplot = setplot

    return claw

    
#--------------------------
    
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

