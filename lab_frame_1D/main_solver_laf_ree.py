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
import reactive_euler_roe
import steadyState as steadyState

plt.close('all')

QOF = [[],[]]
tdata = np.array([])

def b4step(solver, state):
    global QOF, tdata
   # Update shock speed  
    gamma = state.problem_data['gamma']
    qheat = state.problem_data['qheat']    
    Ea = state.problem_data['Ea'] 
    T_ign = state.problem_data['T_ign']  
    gamma1 = gamma-1
    
    rho = state.q[0,:]  
    u = state.q[1,:]/state.q[0,:]
    pressure = gamma1*(state.q[2,:]-0.5*state.q[1,:]**2/state.q[0,:] -
                        qheat*state.q[3,:])
    T = pressure / rho
    animation = state.problem_data['animation'] 
    # Plot
    if animation:
        x=state.grid.x.centers      
        h[0].set_data(x,pressure)
        ax.figure.canvas.draw()
        plt.show()     
    # Compute QOF
    tot_var = np.sum(np.abs(np.diff(rho)))
    QOF[0] = np.append(QOF[0], state.t)
    QOF[1] = np.append(QOF[1], tot_var)
    
#    print(state.q[:,-1])
#    print(state.qbc[:,:])
#    time.sleep(1)

#    print(QOF)    
#    tdata = np.append(tdata, state.t)
#    print(QOF)
#    print(tdata)
#    time.sleep(1)

def total_variation(state):
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
#        print(num_ghost)
#        time.sleep(1)
        gamma = state.problem_data['gamma']
        qheat = state.problem_data['qheat']  
        q = state.q

#        qbc[0,i] = q[0,-1]  
#        qbc[1,i] = q[1,-1]
#        qbc[2,i] = q[2,-1]
#        qbc[3,i] = q[3,-1]        
#        print(qbc[:,:])
        qbc[0,-i-1] = 1  
        qbc[1,-i-1] = 0
        qbc[2,-i-1] = 1/(gamma-1) + qheat
        qbc[3,-i-1] = 1  
#        q = state.q
#        print(i)
#        print(q[:,-1])
#        print(qbc[:,i])
#        time.sleep(2)
        

def main_solver_laf_ree(Ea = 1.65, qheat = 1.7, gamma = 1.4, T_ign = 1.5, 
                        mx=300, xmax=30, tmax=10,dtout=1, 
                        outdir='./_output', iplot=1, animation=0):
                            
        
    gamma1 = gamma-1
    # Set up the solver object
#    solver = pyclaw.ClawSolver1D(reactive_euler_roe)
    solver = pyclaw.ClawSolver1D(reactive_euler_roe)    
    solver.step_source = step_reaction
    solver.before_step = b4step
    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.user_bc_upper = custom_bc
    solver.num_waves = 4

    # Set the domain
    x = pyclaw.Dimension('x',0.0,xmax,mx)
    domain = pyclaw.Domain([x])
    num_eqn = 4

    # Set state    
    state = pyclaw.State(domain,num_eqn)

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
    state.problem_data['fspeed'] = D
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
plt.ylim([0, 2])

gamma = 1.2
gamma1 = gamma-1
q_asympt = 1.7
theta = 2.
qheat = q_asympt*gamma1*gamma
Ea = theta/(gamma1**2)
#qheat = 50
#Ea = 26

tmax = 500
mx = 1000

main_solver_laf_ree(Ea = Ea, qheat = qheat, gamma = gamma, T_ign = 1.05, 
                        mx=mx, xmax=xmax, tmax=tmax, dtout=1, 
                        outdir='./_output', iplot=0, animation=1)   

