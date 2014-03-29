# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 18:47:38 2014

@author: maltezfaria
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 14:11:55 2014

@author: maltezfaria
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 16:27:02 2014

@author: maltezfaria
"""
import numpy as np
import time
import matplotlib.pylab as plt
from clawpack import pyclaw
#import clawpack.petclaw as pyclaw
from clawpack import riemann
import reactive_euler
import steadyState as steadyState

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


QOF = [[],[]]
tdata = np.array([])

def b4step(solver, state):
    global QOF, tdata
    # Update shock speed 
    gamma = state.problem_data['gamma']
    qheat = state.problem_data['qheat'] 
    gamma1 = gamma-1
    q = state.q    
    pressure = gamma1*(q[2,-1]-0.5*q[1,-1]**2/q[0,-1] -
                        qheat*q[3,-1])

    dF1 = q[1,-1] - 0
    dq1 = q[0,-1] - 1
    D1 = dF1/dq1

    dF2 = q[1,-1]**2/q[0,-1] + pressure - 1
    dq2 = q[1,-1]
    D2 = dF2/dq2
    
    dF3 = q[1,-1]*q[2,-1]/q[0,-1] + q[1,-1]/q[0,-1] * pressure - 0
    dq3 = q[2,-1] - (1/(gamma1) + qheat)
    D3 = dF3/dq3
    
    D = (D1+D2+D3)/3    
#    
#    print('D before=' + str(state.problem_data['fspeed']))
    state.problem_data['fspeed'] = D
#    print('D after=' + str(state.problem_data['fspeed']))
#    time.sleep(1)
#    print(D1, D2, D3)
#    print(D)
#    time.sleep(1)
    
    rho = state.q[0,:]        
    # Plot
    animation = state.problem_data['animation'] 
    if animation:
        x=state.grid.x.centers      
        h[0].set_data(x,rho)
        ax.figure.canvas.draw()
        plt.show()     
    # Compute QOF
#    tot_var = np.sum(np.abs(np.diff(rho)))
    QOF[0] = np.append(QOF[0], state.t)
    QOF[1] = np.append(QOF[1], D)

#    print(QOF)    
#    tdata = np.append(tdata, state.t)
#    print(QOF)
#    print(tdata)
#    time.sleep(1)


def total_variation(state):
    import numpy as np

    rho = state.q[0,:]
#    dx=state.grid.delta[0];
    variation = np.abs(np.diff(rho))
    state.F[0,:] = np.append(variation, 0)

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
    
def custom_bc(state,dim,t,qbc,num_ghost):
    for i in xrange(num_ghost):
#        print(i)
        time.sleep(1)
        gamma = state.problem_data['gamma']
        qheat = state.problem_data['qheat']
        gamma1 = gamma-1
        q = state.q
        
        qbc[0,-i-1] = q[0,-1]  
        qbc[1,-i-1] = q[1,-1]
        qbc[2,-i-1] = q[2,-1]
        qbc[3,-i-1] = 1

def main_solver_saf_ree(Ea = 1.65, qheat = 1.7, gamma = 1.4, T_ign = 1.1, 
                        mx=300, xmax=30, tmax=10,dtout=1, 
                        outdir='./_output', iplot=1, animation=0):
                            
        
    gamma1 = gamma-1
    # Set up the solver object
    solver = pyclaw.ClawSolver1D(reactive_euler)
    solver.step_source = step_reaction
    solver.before_step = b4step
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.custom
    solver.user_bc_upper = custom_bc
    


    
    # Set the domain
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 4

    # Set state    
    state = pyclaw.State(domain,num_eqn)


    # Set initial condition
    x =state.grid.x.centers
#    xs = xmax-5
    xs = np.inf    
    xdet = x[x<xs]
    rhol, Ul, pl, laml, D0, k = steadyState.steadyState(qheat, Ea, gamma, xdet)
    ul = Ul + D0
    Yl = 1-laml
    rhor = 1.*np.ones(np.shape(x[x>=xs]))
    ur = 0.0*np.ones(np.shape(x[x>=xs]))
    pr = 1*np.ones(np.shape(x[x>=xs]))
    Yr = 1.*np.ones(np.shape(x[x>=xs]))
    rho = np.append(rhol, rhor)
    u = np.append(ul, ur)
    p = np.append(pl, pr)
    Y = np.append(Yl, Yr)
    
    state.q[0,:] = rho
    state.q[1,:] = rho * u
    state.q[2,:] = p/gamma1 + rho*u**2/2 + qheat * rho * Y
    state.q[3,:] = rho * Y
    state.mF = 1

    # Set initial condition
#    x =state.grid.x.centers
#    rho, U, p, lam, D, k = steadyState.steadyState(qheat, Ea, gamma, x)
#    u = U + D
#    Y = 1-lam
#    
#    state.q[0,:] = rho
#    state.q[1,:] = rho * u
#    state.q[2,:] = p/gamma1 + rho*u**2/2 + qheat * rho * Y
#    state.q[3,:] = rho * Y
#    state.mF = 1
#

    # Fill dictory of problem data
    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['qheat']= qheat
    state.problem_data['Ea'] = Ea
    state.problem_data['T_ign'] = T_ign
    state.problem_data['fspeed'] = D0
    state.problem_data['k'] = k
    state.problem_data['animation'] = animation 




#    plt.plot(x, u)
#    plt.show()

    # Set up controller   

    claw = pyclaw.Controller()
    claw.compute_F = total_variation
    claw.tfinal = tmax
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    #Set up data writting and plotting 
    claw.output_style = 1
    nout = np.ceil(tmax/dtout)
    claw.num_output_times = nout    
    claw.setplot = setplot
    claw.outdir = outdir
    claw.keep_copy = iplot    
    
    #Run the simulation
    claw.run()    

    claw.plot()
    
    #Save QOF
    np.save('qof',QOF)
    
    return QOF

    
xmax=50
plt.ion()
ax = plt.gca()    
h = ax.plot([],[],'ro-')
plt.xlim([0, xmax])
plt.ylim([0, 3])

main_solver_saf_ree(Ea = 4, qheat = 1, gamma = 1.2, T_ign = 1.1, 
                        mx=1000, xmax=xmax, tmax=100, dtout=1, 
                        outdir='./_output', iplot=1, animation=1)   
