#!/usr/bin/env python
# encoding: utf-8
"""
Shock tube test problem
===================================

Solve the one-dimensional Euler equations for inviscid, compressible flow:

.. math::
    \rho_t + (\rho u)_x & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x & = 0 \\
    E_t + (u (E + p) )_x & = 0.

The fluid is an ideal gas, with pressure given
 by :math:`p=\rho (\gamma-1)e` where
e is internal energy.
"""
gamma = 1.4
gamma1 = gamma - 1.


def setup(use_petsc=False, outdir='./_output', solver_type='classic'):

    from clawpack import pyclaw
    from clawpack import riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann.euler_with_efix_1D)
        solver.weno_order = 5

        #solver.lim_type = 2
        solver.char_decomp = 1
        solver.time_integrator = 'SSP104'
    else:
        solver = pyclaw.ClawSolver1D(riemann.euler_with_efix_1D)

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.limiters = 1


    # Initialize domain
    mx=200;
    x = pyclaw.Dimension('x',0.0,2.0,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 3
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1


    rhol =  1.
    rhor = 0.125
    ul = 0.
    ur = 0.
    pl = 1.
    pr = 0.1
    xs1 = 0.75    
    xs2 = 1.25     
    
    x =state.grid.x.centers
    state.q[0,:] = (x>xs1)*(x<xs2)*rhol + ~((x>xs1)*(x<xs2))*rhor
    state.q[1,:] = (x>xs1)*(x<xs2)*rhol*ul + ~((x>xs1)*(x<xs2))*rhor*ur
    state.q[2,:] = ((x>xs1)*(x<xs2)*(pl/gamma1 + rhol*ul**2/2) + 
                    ~((x>xs1)*(x<xs2))*(pr/gamma1 + rhor*ur**2/2) )

#    rhol =  1.
#    rhor = 0.125
#    ul = 0.
#    ur = 0.
#    pl = 1.
#    pr = 0.1
#    xs = 0.5
#
#    x =state.grid.x.centers
#    state.q[0,:] = (x<xs)*rhol + (x>=xs)*rhor
#    state.q[1,:] = (x<xs)*rhol*ul + (x>=xs)*rhor*ur
#    state.q[2,:] = (x<xs)*(pl/gamma1 + rhol*ul**2/2) + (x>=xs)*(pr/gamma1 + rhor*ur**2/2)
#
#
    claw = pyclaw.Controller()
    claw.tfinal = 0.3
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


def pressure(current_data):
    """Computes the pressure from the conserved quantities"""
    pressure = gamma1*(current_data.q[2,:]-
                0.5*current_data.q[1,:]**2/current_data.q[0,:])
    return pressure

def temperature(current_data):
    return pressure(current_data)/current_data.q[0,:]

def velocity(current_data):
    """Computes the velocity from the conserved quantities"""
    velocity = current_data.q[1,:]/current_data.q[0,:]
    return velocity

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
