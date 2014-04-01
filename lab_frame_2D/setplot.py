def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    from clawpack.visclaw import colormaps


    def label_axes(current_data):
        import matplotlib.pyplot as plt
        plt.xlabel('x')
        plt.ylabel('y')

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Pressure plot
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Density'
    plotaxes.scaled = False      # so aspect ratio is 1
    plotaxes.afteraxes = label_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem.plot_var = 0
    plotitem.add_colorbar = False
    

    # Tracer plot
    plotfigure = plotdata.new_plotfigure(name='Tracer', figno=1)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Tracer'
    plotaxes.scaled = False     # so aspect ratio is 1
    plotaxes.afteraxes = label_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax=1.0
    plotitem.plot_var = 4
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = False
    

    # Energy plot
    plotfigure = plotdata.new_plotfigure(name='Energy', figno=2)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'V'
    plotaxes.scaled = False      # so aspect ratio is 1
    plotaxes.afteraxes = label_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
#    plotitem.pcolor_cmin = 2.
#    plotitem.pcolor_cmax=18.0
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    
    return plotdata
