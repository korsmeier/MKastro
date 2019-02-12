import  matplotlib                  as      mpl
mpl.use('Agg')
import  matplotlib.pyplot           as      plt
import  matplotlib.patches          as      patches
import  matplotlib.colors           as      colors

from    matplotlib.path             import  Path
from    matplotlib.colors           import  colorConverter
from    matplotlib.legend_handler   import  HandlerLine2D
from    matplotlib                  import  rc
import  matplotlib.ticker           as      tck
from    matplotlib.ticker           import  FuncFormatter


def new_plot( xlabel, ylabel, xscale='linear', yscale='linear', print_size=10, label_size=25, sizex=1.0, sizey=0.75):
    
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    
    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)
    
    plot.set_xlabel ( xlabel )
    plot.set_ylabel ( ylabel )
    
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )
    
    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

    
    return plot, fig



def new_plot_res( xlabel='', ylabel='', xscale='log', yscale='log', print_size=10, label_size=25, sizex=1.0, sizey=1.0):
    
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    
    fig          = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot_main    = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    plot_res     = plt.subplot2grid((4, 1), (3, 0), rowspan=1, sharex=plot_main)
    
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)

    
    plot_res .set_xlabel ( xlabel )
    plot_main.set_ylabel ( ylabel )
    
    plot_main.set_xscale ( xscale   )
    plot_main.set_yscale ( yscale   )
    plot_res .set_yscale ( 'linear' )
    
    plot_main.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot_main.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    
    plot_res. tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot_res. tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    
    plot_main.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    plot_res .grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    
    plt.setp( plot_main.get_xticklabels(), visible=False )
    
    return plot_main, plot_res, fig

