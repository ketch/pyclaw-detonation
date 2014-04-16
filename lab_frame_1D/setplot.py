def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    from clawpack.visclaw import colormaps
    import numpy as np
    from znd_wave_2D import gamma

    #Define pressure for plotting
    # def pressure(current_data):
    #     qheat  = 5
    #     q = current_data.q
    #     rho = q[0,:,:]
    #     u   = q[1,:,:]/rho
    #     v   = q[2,:,:]/rho
    #     press = gamma1 *  (q[3,:,:] - 0.5*rho*(u**2 + v**2) - qheat*q[4,:,:])
    #     return press
        
    # def temperature(current_data):
    #     rho = current_data.q[0,:,:]
    #     press = pressure(current_data)
    #     temper = np.sqrt(press/rho)
    #     return temper

    def y_velocity(current_data):
        return current_data.q[2,:,:] / current_data.q[0,:,:]

    def reacVar(current_data):
        #print(vars(current_data))
        #print(current_data.plotdata)
        return current_data.q[4,:,:] / current_data.q[0,:,:]
        
    def label_axes(current_data):
        import matplotlib.pyplot as plt
        plt.xlabel('x')
        plt.ylabel('y')
        #plt.xlim((20,28))

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Density plot
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Density'
    plotaxes.scaled = False      # so aspect ratio is 1
    plotaxes.afteraxes = label_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem.plot_var = 0
    plotitem.add_colorbar = True
    

    # Tracer plot
    plotfigure = plotdata.new_plotfigure(name='Tracer', figno=1)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Tracer'
    plotaxes.scaled = False     # so aspect ratio is 1
    plotaxes.afteraxes = label_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 1.0
    plotitem.plot_var = reacVar
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    

    # y velocity
    plotfigure = plotdata.new_plotfigure(name='V', figno=2)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'V'
    plotaxes.scaled = False      # so aspect ratio is 1
    plotaxes.afteraxes = label_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
#    plotitem.pcolor_cmin = 2.
#    plotitem.pcolor_cmax=18.0
    plotitem.plot_var = y_velocity
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    
    return plotdata
