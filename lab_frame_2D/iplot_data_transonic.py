from clawpack import petclaw
from setplot_transonic import setplot_transonic

#ip = petclaw.plot.interactive_plot(setplot=setplot,outdir="_output_weak_heat_release1")
ip = petclaw.plot.interactive_plot(setplot=setplot_transonic,outdir="_output_transonic")
