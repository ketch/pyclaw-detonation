#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import mpi4py.MPI as mpi
from setplot import setplot
#from pprint import pprint

gamma = 1.2
gamma1 = gamma - 1.
q_asympt = 1.7
theta = 1.65

gamma1a = gamma;
epsi = (gamma-1)/(gamma1a);
Ea = theta/epsi**2*gamma/gamma1a;
qheat = epsi*q_asympt;

qheat = 0.4
Ea = 50.

T_ign = 0

xmax = 10.
ymax = 1
# ymax = 20
mx = 400
my = ymax/xmax*mx
tfinal = 4
num_output_times=10

def b4step(solver, state):
    pass

def qinit(state,domain):

    grid = state.grid
    ((ixlower,ixupper),(iylower,iyupper)) = domain.patch._da.getRanges()
    #print(ixlower, ixupper)

    dx = grid.delta[0] #Could be problem for non-uniform grid
    xglobal = np.linspace(dx/2,xmax-dx/2,mx)

    x =grid.x.centers
    y =grid.y.centers


    rhol =  1.
    rhor = 0.125
    ul = 0.5
    ur = 0.5
    pl = 1.
    pr = 0.1
    xs1 = -1.
    xs2 = 3.

    for j in range(y.size):
        state.q[0,:,j] = (x>xs1)*(x<xs2)*rhol + ~((x>xs1)*(x<xs2))*rhor
        state.q[1,:,j] = (x>xs1)*(x<xs2)*rhol*ul + ~((x>xs1)*(x<xs2))*rhor*ur
        state.q[2,:,j] = 0
        state.q[3,:,j] = ((x>xs1)*(x<xs2)*(pl/gamma1 + rhol*ul**2/2) +
                    ~((x>xs1)*(x<xs2))*(pr/gamma1 + rhor*ur**2/2) )
        state.q[4,:,j] = 1

    state.problem_data['xfspeed']= 0
    state.problem_data['k']= 0

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
          outdir='_output_transonic', disable_output=False, mx=mx, my=my, tfinal=tfinal,
          num_output_times = num_output_times):
    """
    Solve the reactive Euler equations of compressible fluid dynamics.
    """
    from clawpack import riemann

    import reactive_euler_roe_2D
    import reactive_euler_efix_roe_2D

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw


    solver = pyclaw.ClawSolver2D(reactive_euler_efix_roe_2D)
    #solver = pyclaw.ClawSolver2D(riemann.euler_5wave_2D)
    solver.order = 1
    solver.limiters = [4,4,4,4,2]
    solver.step_source=step_Euler_reaction
    solver.before_step = b4step

    # Initialize domain
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    y = pyclaw.Dimension('y',0.0,ymax,my)
    domain = pyclaw.Domain([x,y])
    solver.num_eqn = 5
    solver.num_waves = 5

    state = pyclaw.State(domain,solver.num_eqn)
    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat

    # The k and xfspeed problem data are set inside qinit
    qinit(state,domain)

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.dt_initial=0.005
    solver.source_split = 1
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.wall
    solver.bc_upper[1]=pyclaw.BC.wall

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    #claw.output_format = ['hdf5', 'ascii']
    claw.keep_copy = False

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
