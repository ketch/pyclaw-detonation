def setplot(plotdata):
#--------------------------
    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    from clawpack.visclaw import colormaps
    import numpy as np
    from znd_wave_1D import gamma, gamma1, qheat
    
    def pressure(current_data):
        q = current_data.q
        rho = q[0,:]
        u   = q[1,:]/rho
        press = gamma1 *  (q[2,:] - 0.5*rho*(u**2) - qheat*q[3,:])
        return press

    def temperature(current_data):
        rho = current_data.q[0,:]
        press = pressure(current_data)
        temper = np.sqrt(press/rho)
        return temper

    def reacVar(current_data):
        return current_data.q[3,:] / current_data.q[0,:]
     

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'pressure'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = pressure
    plotitem.plotstyle = '-'
    plotitem.color = 'b'


    # Figure for Y
    plotfigure = plotdata.new_plotfigure(name='Y', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Y'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = reacVar
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

    # Figure for temperature
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
        
