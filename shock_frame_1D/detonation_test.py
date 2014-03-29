#!/usr/bin/env python
# encoding: utf-8
"""
Shock tube test problem
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
gamma = 1.4
gamma1 = gamma - 1.
qheat = 1.
Ea = 1.
k = -1
T_ign = 0.8
fspeed = 0

def setup(outdir='./_output'):

    import numpy as np
    from clawpack import pyclaw
    from clawpack import riemann
    import reactive_euler 

    solver = pyclaw.ClawSolver1D(reactive_euler)
    solver.step_source = step_reaction

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap

    # Initialize domain
    mx=2000;
    x = pyclaw.Dimension('x',0.0,200.0,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 4
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat
    state.problem_data['Ea'] = Ea
    state.problem_data['k'] = k
    state.problem_data['T_ign'] = T_ign
    state.problem_data['fspeed'] = fspeed

    rhol =  1.
    rhor = 0.125
    ul = 0.5
    ur = 0.
    pl = 1.
    pr = 0.1
    Yl = 0
    Yr = 1
    xs = 20
    
    x =state.grid.x.centers
    state.q[0,:] = (x<xs)*rhol + (x>=xs)*rhor
    state.q[1,:] = (x<xs)*rhol*ul + (x>=xs)*rhor*ur
    state.q[2,:] = (x<xs)*(pl/gamma1 + rhol*ul**2/2) + (x>=xs)*(pr/gamma1 + rhor*ur**2/2)
    state.q[3,:] = (x<xs)*(rhol*Yl) + (x>=xs)*(rhor*Yr)


    claw = pyclaw.Controller()
    claw.tfinal = 100
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    #state.mp = 1
    #claw.compute_p = pressure
    #claw.write_aux_always = 'True'

    return claw

def step_reaction(solver, state, dt):
    q = state.q
    q[0,:] = q[0,:]
    q[1,:] = q[1,:]
    q[2,:] = q[2,:]
    q[3,:] = q[3,:] + dt*omega(q)


def omega(q):
    """ Reaction rate function
    """
    import numpy as np
    pressure = gamma1*(q[2,:]-0.5*q[1,:]**2/q[0,:] - 
                        qheat*q[3,:])
    T = pressure / q[0,:]
    return k*(q[3,:])*np.exp(-Ea / T )*(T>T_ign)

def temperature(current_data):
    pressure = gamma1*(current_data.q[2,:]-0.5*current_data.q[1,:]**2/current_data.q[0,:] - 
                        qheat*current_data.q[3,:])
    T = pressure / current_data.q[0,:]
    return T

def pressure(current_data):
    """Computes the pressure from the conserved quantities"""
    pressure = gamma1*(current_data.q[2,:]-0.5*current_data.q[1,:]**2/current_data.q[0,:] - 
                        qheat*current_data.q[3,:])
    return pressure

def velocity(current_data):
    """Computes the velocity from the conserved quantities"""
    velocity = current_data.q[1,:]/current_data.q[0,:]
    return velocity

def reacVar(current_data):
    """Computes reaction progress variable lambda"""
    out = current_data.q[3,:]/current_data.q[0,:]
    return out

def add_true_solution(current_data):
    '''Adds a plot of true solution
    '''
    from pylab import plot
    x = current_data.x
    t = current_data.t
    ptrue = x*t
    plot(x,ptrue,'r')

#--------------------------
def setplot(plotdata):
#--------------------------
    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'pressure'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = pressure
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

#    plotaxes.afteraxes = add_true_solution

    # Figure for velocity
    plotfigure = plotdata.new_plotfigure(name='Velocity', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = velocity
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

    # Figure for Y
    plotfigure = plotdata.new_plotfigure(name='Y', figno=3)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Y'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = reacVar
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

    # Figure for q[1]
    plotfigure = plotdata.new_plotfigure(name='Temperature', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Temperature'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = temperature
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
