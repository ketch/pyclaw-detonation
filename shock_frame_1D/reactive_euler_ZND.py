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
gamma = 1.2
gamma1 = gamma - 1.
qheat = 5.
Ea = 10
T_ign = 1.1


def compute_f(state):
    import numpy as np
#    print('Entered compute_f')
#    print(state.q[0,:])
    state.F[0] = state.q[0,:]
    
def total_variation(state):
    import numpy as np

    rho = state.q[0,:]
#    dx=state.grid.delta[0];
    variation = np.abs(np.diff(rho))
    state.F[0,:] = np.append(variation, 0)

def step_reaction(solver, state, dt):
    global k
    q = state.q
    q[0,:] = q[0,:]
    q[1,:] = q[1,:]
    q[2,:] = q[2,:]
    q[3,:] = q[3,:] + dt*omega(q, k)


def omega(q, k):
    """ Reaction rate function
    """
    import numpy as np
    pressure = gamma1*(q[2,:]-0.5*q[1,:]**2/q[0,:] -
                        qheat*q[3,:])
    T = pressure / q[0,:]
    return -k*(q[3,:])*np.exp(Ea*(1-1/T))*(T>T_ign)

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


def setup(outdir='./_output'):

    global k

    import numpy as np
    import matplotlib.pylab as plt
    from clawpack import pyclaw
    from clawpack import riemann
    import reactive_euler
    import steadyState as steadyState

    solver = pyclaw.ClawSolver1D(reactive_euler)
    solver.step_source = step_reaction

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap

    # Initialize domain
    mx=2000
    xmax=50
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 4
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat
    state.problem_data['Ea'] = Ea
    state.problem_data['T_ign'] = T_ign

    x =state.grid.x.centers
#    xs = xmax-5
    xs = np.inf    
    xdet = x[x<xs]

    rhol, Ul, pl, laml, D, k = steadyState.steadyState(qheat, Ea, gamma, xdet)
    ul = Ul + D
    Yl = 1-laml
    state.problem_data['fspeed'] = D
#    state.problem_data['k'] = xscale

    rhor = 1.*np.ones(np.shape(x[x>=xs]))
    ur = 0.0*np.ones(np.shape(x[x>=xs]))
    pr = 1*np.ones(np.shape(x[x>=xs]))
    Yr = 1.*np.ones(np.shape(x[x>=xs]))

    rho = np.append(rhol, rhor)
    u = np.append(ul, ur)
    p = np.append(pl, pr)
    Y = np.append(Yl, Yr)

#    plt.plot(x, rho, x, Y, x, u, x, p)
#    plt.legend(('rho', 'Y', 'u','p'))
#    plt.show()

    state.q[0,:] = rho
    state.q[1,:] = rho * u
    state.q[2,:] = p/gamma1 + rho*u**2/2 + qheat * rho * Y
    state.q[3,:] = rho * Y

    state.mF = 1

    claw = pyclaw.Controller()
    claw.compute_F = total_variation
    claw.tfinal = 50
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
