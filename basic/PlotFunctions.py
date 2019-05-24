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




#colors      = ['#C0392B','#E74C3C','#9B59B6','#8E44AD','#2980B9','#3498DB','#1ABC9C','#16A085','#27AE60','#2ECC71','#F1C40F','#F39C12','#D35400','#95A5A6','#7F8C8D','#34495E','#2C3E50']
colors      = ['#C0392B','#E74C3C','#9B59B6','#8E44AD','#2980B9','#3498DB','#1ABC9C','#16A085','#196F3D','#1D8348','#F1C40F','#F39C12','#D35400','#95A5A6','#7F8C8D','#34495E','#2C3E50']
colors_dark = ['#7B241C','#943126','#633974','#5B2C6F','#1A5276','#21618C','#117864','#0E6655','#196F3D','#1D8348','#9A7D0A','#9C640C','#873600','#5F6A6A','#515A5A','#212F3C','#1C2833']

cmap_rainbow = mpl.cm.get_cmap('jet')

dashes      = [  (), (2,2), (5,3), (7,2,2,2)   ]

markers     = [  'o','v','^', 's', 'd', '<', '>', 'P', 'X' ]

