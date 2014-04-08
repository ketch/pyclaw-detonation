#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import mpi4py.MPI as mpi
from setplot import setplot
#from pprint import pprint

gamma = 1.2
gamma1 = gamma - 1.
q_asympt = 1.7
theta = 1.8
qheat = q_asympt*gamma1*gamma
Ea = theta/(gamma1**2)

#qheat = 50
#Ea = 20

if mpi.COMM_WORLD.Get_rank()==0:
    print('Heat release Q = {} and activation energy = {}'.format(qheat,Ea)) 

T_ign = 1.005

xmax = 30.
ymax = 50.
mx = 200
my = 400
xs = 25.
tfinal = 1000.
num_output_times=500

def b4step(solver, state):
    pass

def qinit(state,domain,xs=xs):
    import steadyState as steadyState

    grid = state.grid
    ((ixlower,ixupper),(iylower,iyupper)) = domain.patch._da.getRanges()
    #print(ixlower, ixupper)

    dx = grid.delta[0] #Could be problem for non-uniform grid
    xglobal = np.linspace(dx/2,xmax-dx/2,mx)

    x =grid.x.centers
    y =grid.y.centers
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

    for i in range(y.size):
        #pprint(y.size)
        pert = 0.0001 * np.sin(2*4*np.pi*y[i]/ymax) * (x<xs)
        #pert=0
        #print(np.shape(rho),np.shape(x))
        state.q[0,:,i] = rho + pert
        state.q[1,:,i] = rho * u +pert
        state.q[2,:,i] = 0. + pert
        state.q[3,:,i] = p/gamma1 + rho*(u**2 + v**2)/2 + qheat * rho * Y  + pert
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

def custom_bc(state,dime,t,qbc,num_ghost):
    """Right boundary condition
    """
    rho_a = 1
    p_a = 1
    u_a = 0
    Y_a = 1
    
    

def setup(use_petsc=False,kernel_language='Fortran',solver_type='classic',
          outdir='_output', disable_output=False, mx=mx, my=my, tfinal=tfinal,
          num_output_times = num_output_times):
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
    qinit(state,domain)

    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.dt_initial=0.005
    solver.source_split = 1
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.custom
    solver.user_bc_upper= custom_bc
    solver.bc_lower[1]=pyclaw.BC.periodic
    solver.bc_upper[1]=pyclaw.BC.periodic

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


