import  math
import  os
import  numpy                   as      np
from    operator                import  itemgetter

import  matplotlib              as      mpl
mpl.use('Agg')
import  matplotlib.pyplot       as      plt
import  matplotlib.patches      as      mpatches
from    matplotlib.colors       import  BoundaryNorm
from    matplotlib.colors       import  colorConverter
from    matplotlib.ticker       import  MaxNLocator
from    matplotlib.path         import  Path
from    matplotlib              import  rc

fontset='dejavuserif'

def prepare_triangle(npar, label_size=None, print_size=15, **kwargs):
    #
    if label_size==None:
        label_size = int(2*print_size/np.sqrt(npar))
    #
    plt.close('all')
    #
    # General settings: fontsize and plotgrid
    #
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 1
    mpl.rcParams['mathtext.fontset']=fontset
    #
    if not 'max_n_locator' in kwargs:
        loc = 6
        if npar > 3:
            loc = 5
        if npar > 5:
            loc = 4
        if npar > 7:
            loc = 3
        kwargs['max_n_locator'] = loc
    #
    dim = npar
    fig, plot = plt.subplots(dim, dim)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    fig.set_size_inches(print_size, print_size)

    if dim==1:
        plotArray=[]
        plotArray.append([plot])
    else:
        plotArray=plot

    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop over all triangle:  remove all axis  #
    # * * * * * * * * * * * * * * * * * * * * * #
    for i in range( 0, len( plotArray )  ):
        for j in range( 0, len( plotArray[i] )  ):
            plotArray[i][j].axis('off')

    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop triangles:                           #
    #    -set labels and ticks                  #
    # * * * * * * * * * * * * * * * * * * * * * #
    for i in range(0, npar):
        for j in range(0, i+1):
            iP = i #-1
            jP = j
            plotArray[iP][jP].axis('on')
            #
            for tick in plotArray[iP][jP].get_xticklabels():
                tick.set_rotation(90)
                tick.set_horizontalalignment('right')
            plotArray[iP][jP].xaxis.set_major_locator(MaxNLocator(kwargs['max_n_locator']))
            plotArray[iP][jP].yaxis.set_major_locator(MaxNLocator(kwargs['max_n_locator']))
            plotArray[iP][jP].xaxis.set_tick_params(length=6, width=1, top=True, right=True, direction='in', pad=10)
            plotArray[iP][jP].yaxis.set_tick_params(length=6, width=1, top=True, right=True, direction='in', pad=10)
            if j==0:
                plotArray[iP][jP].set_ylabel( 'par %i' %i )
            if i==npar-1:
                plotArray[iP][jP].set_xlabel( 'par %i' %j )
            if i<npar-1:
                plotArray[iP][jP].xaxis.set_tick_params( labelsize=0, top=True, right=True, direction='in', pad=10  )
                plotArray[iP][jP].xaxis.get_offset_text().set_fontsize(0)
            if j>0:
                plotArray[iP][jP].yaxis.set_tick_params( labelsize=0, top=True, right=True, direction='in', pad=10  )
                plotArray[iP][jP].yaxis.get_offset_text().set_fontsize(0)
            #
    return fig, plotArray

def set_triangle_labels_and_ranges( plotArray, labels, ranges=[] ):
    npar = len(plotArray)
    y_offset = 0
    if len(labels)==npar+1:
        y_offset = +1
    for i in range(0, npar):
        for j in range(0, i+1):
            iP = i 
            jP = j
            if jP==0:
                plotArray[iP][jP].set_ylabel( labels[i+y_offset] )
            if iP==npar-1:
                plotArray[iP][jP].set_xlabel( labels[j] )
            if len(ranges):
                plotArray[iP][jP].set_xlim( ranges[j] )
                if iP!=jP:
                    plotArray[iP][jP].set_ylim( ranges[i+y_offset] )
    if len(labels)==npar:
        plotArray[0][0].set_ylabel('$-2\log(\mathcal{L})$')

def unify_ranges_from_diagonal( plotArray ):
    npar = len(plotArray)
    ranges = []
    for i in range(0, npar):
        ranges.append( plotArray[i,i].get_xlim() )
    for i in range(0, npar):
        for j in range(0, i+1):
            iP = i
            jP = j
            if len(ranges):
                plotArray[iP][jP].set_xlim( ranges[j] )
                if iP!=jP:
                    plotArray[iP][jP].set_ylim( ranges[i] )


def set_triangle_scales( fig, plotArray, scales, **kwargs ):
    npar = len(plotArray)
    y_offset = 0
    aax2 = []
    for i in range(0, npar):
        l = []
        for j in range(0, npar):
            l.append(  fig.add_axes(plotArray[i][j].get_position())  )
            if j>i:
                l[-1].axis('off')
        aax2.append(l)
    if not 'max_n_locator' in kwargs:
        loc = 6
        if npar > 3:
            loc = 5
        if npar > 5:
            loc = 4
        if npar > 7:
            loc = 3
        kwargs['max_n_locator'] = loc
    if len(scales)==npar+1:
        y_offset = +1
    for i in range(0, npar):
        for j in range(0, i+1):
            iP = i
            jP = j
            plotArray[iP][jP].axis('off')
            aax2[iP][jP].patch.set_alpha(0.0)
            aax2[iP][jP].set_xlim( plotArray[iP][jP].get_xlim() )
            aax2[iP][jP].set_ylim( plotArray[iP][jP].get_ylim() )
            if scales[j]=='log':
                aax2[iP][jP].set_xlim( np.power(10, np.array(plotArray[iP][jP].get_xlim())) )
            aax2[iP][jP].set_xscale( scales[j] )
            if iP!=jP:
                if scales[i]=='log':
                    aax2[iP][jP].set_ylim( np.power(10, np.array(plotArray[iP][jP].get_ylim())) )
                aax2[iP][jP].set_yscale( scales[i+y_offset] )
            aax2[iP][jP].xaxis.set_tick_params(length=6, width=1, top=True, right=True, direction='in', pad=10)
            aax2[iP][jP].yaxis.set_tick_params(length=6, width=1, top=True, right=True, direction='in', pad=10)
            if scales[j]=='log':
                aax2[iP][jP].xaxis.set_tick_params('minor', length=3, width=1, top=True, right=True, direction='in', pad=10)
            if iP!=jP:
                if scales[i]=='log':
                    aax2[iP][jP].yaxis.set_tick_params('minor', length=3, width=1, top=True, right=True, direction='in', pad=10)
            for tick in aax2[iP][jP].get_xticklabels(which='both'):
                tick.set_rotation(90)
            if j==0:
                aax2[iP][jP].set_ylabel( plotArray[iP][jP].get_ylabel() )
            if i==npar-1:
                aax2[iP][jP].set_xlabel( plotArray[iP][jP].get_xlabel() )
            if i<npar-1:
                aax2[iP][jP].xaxis.set_tick_params( 'both', labelsize=0, top=True, right=True, direction='in', pad=10  )
                aax2[iP][jP].xaxis.get_offset_text().set_fontsize(0)
            if j>0:
                aax2[iP][jP].yaxis.set_tick_params( 'both', labelsize=0, top=True, right=True, direction='in', pad=10  )
                aax2[iP][jP].yaxis.get_offset_text().set_fontsize(0)
            if len(scales)==npar and iP==0 and jP==0:
                aax2[iP][jP].set_ylabel('$-2\log(\mathcal{L})$')
    return aax2


def draw_triangle( fig, plotArray, draw_function, vector_chi2, matrix_parameter, vector_positions, **kwargs_draw_function ):
    im = None
    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop triangles:                           #
    #    - plot in lower triangle               #
    # * * * * * * * * * * * * * * * * * * * * * #
    if( len(vector_positions)>len(plotArray[:,0]) ):
        y_offset=-1
    else:
        y_offset=0
    x_offset=0
    for i in range(0, len(vector_positions)):
        if vector_positions[i]<0:
            continue
        for j in range(0, i):
            if vector_positions[j]<0:
                continue
            #
            iP = i + y_offset#-1
            jP = j + x_offset
            #
            plotArray[iP][jP].axis('on')
            im = draw_function(fig, plotArray[iP][jP], vector_chi2, matrix_parameter[:,[vector_positions[j],vector_positions[i]]], vector_positions=[vector_positions[j],vector_positions[i]], **kwargs_draw_function)
    return im


def draw_diagonal( fig, plotArray, draw_function, vector_chi2, matrix_parameter, vector_positions, **kwargs_draw_function ):
    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop diagonal :                           #
    #    - plot on diagonal                     #
    # * * * * * * * * * * * * * * * * * * * * * #
    for i in range(0, len(vector_positions)):
        if vector_positions[i]<0:
            continue
        #
        iP = i #-1
        jP = i
        plotArray[iP][jP].axis('on')
        draw_function(fig, plotArray[iP][jP], vector_chi2, matrix_parameter[:,vector_positions[i]], **kwargs_draw_function)


from   numpy import linalg as LA
# 1 sigma contour of a cov matrix V
def contour_deta_chi2_eq_1(V,i,j):
    M = LA.inv(V[i:j+1:j-i,i:j+1:j-i])
    def x_vec(t):
        w,v = LA.eigh(M)
        return np.sqrt(1./w[0])*(v[0])[0]*np.cos(t)+np.sqrt(1./w[1])*(v[1])[0]*np.sin(t), np.sqrt(1./w[0])*(v[0])[1]*np.cos(t)+np.sqrt(1./w[1])*(v[1])[1]*np.sin(t)
    t = np.arange(0,2*math.pi+0.11,0.1)
    return x_vec(t)

def draw_triangle_covariance_matrix( fig, plotArray, mean, V, vector_positions, color='blue', dashes=(), lw=1.0 ):
    im = None
    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop triangles:                           #
    #    - plot in lower triangle               #
    # * * * * * * * * * * * * * * * * * * * * * #
    if( len(vector_positions)>len(plotArray[:,0]) ):
        y_offset=-1
    else:
        y_offset=0
    x_offset=0
    for i in range(0, len(vector_positions)):
        if vector_positions[i]<0:
            continue
        for j in range(0, i):
            if vector_positions[j]<0:
                continue
            #
            iP = i + y_offset#-1
            jP = j + x_offset
            #
            plotArray[iP][jP].axis('on')
            line = contour_deta_chi2_eq_1(V,vector_positions[j],vector_positions[i])
        
            plotArray[iP][jP].plot( line[0] * np.sqrt(2.3) + mean[j] , line[1]  * np.sqrt(2.3) + mean[i], color=color, dashes=dashes, lw=lw )
    return im
