import  numpy                   as      np
import  math

import  matplotlib              as      mpl
mpl.use('Agg')
import  matplotlib.pyplot       as      plt
import  matplotlib.patches      as      mpatches
from    matplotlib.patches      import  Rectangle


def get_bayes_in_adaptive_grid( vec_x, vec_y, vec_weight, range_x=None, range_y=None, r=0, r_max=6, r_min=3, min_Npoints=15):
    if range_x==None:
        range_x = [np.amin(vec_x), np.amax(vec_x) ]
    if range_y==None:
        range_y = [np.amin(vec_y), np.amax(vec_y) ]
    x1 = range_x[0]
    x2 = 0.5 *(range_x[0]+range_x[1])
    x3 = range_x[1]
    #
    y1 = range_y[0]
    y2 = 0.5 *(range_y[0]+range_y[1])
    y3 = range_y[1]
    if not r<r_min:
        if (r>=r_max or len(vec_x)<min_Npoints):
            rel_size = np.power(4.,r_max-r)
            if len(vec_weight):
                return [[ x1, x3, y1, y3, (np.array(vec_weight)).sum()/rel_size, rel_size]]
            return [[ x1, x3, y1, y3, 0, rel_size]]
    #
    ur_vx   = []
    ur_vy   = []
    ur_chi2 = []
    #
    ul_vx   = []
    ul_vy   = []
    ul_chi2 = []
    #
    lr_vx   = []
    lr_vy   = []
    lr_chi2 = []
    #
    ll_vx   = []
    ll_vy   = []
    ll_chi2 = []
    #
    for i,x in enumerate(vec_x):
        y    = vec_y   [i]
        weigth = vec_weight[i]
        if x < x1 or x > x3 or y < y1 or y>y3:
            continue
        if x < x2:
            if y < y2:
                ll_vx.append(x)
                ll_vy.append(y)
                ll_chi2.append(weigth)
            else:
                ul_vx.append(x)
                ul_vy.append(y)
                ul_chi2.append(weigth)
        else:
            if y < y2:
                lr_vx.append(x)
                lr_vy.append(y)
                lr_chi2.append(weigth)
            else:
                ur_vx.append(x)
                ur_vy.append(y)
                ur_chi2.append(weigth)
    #
    l1 = get_bayes_in_adaptive_grid( ll_vx, ll_vy, ll_chi2, [x1,x2], [y1,y2], r+1, r_max, r_min, min_Npoints )
    l2 = get_bayes_in_adaptive_grid( ul_vx, ul_vy, ul_chi2, [x1,x2], [y2,y3], r+1, r_max, r_min, min_Npoints )
    l3 = get_bayes_in_adaptive_grid( lr_vx, lr_vy, lr_chi2, [x2,x3], [y1,y2], r+1, r_max, r_min, min_Npoints )
    l4 = get_bayes_in_adaptive_grid( ur_vx, ur_vy, ur_chi2, [x2,x3], [y2,y3], r+1, r_max, r_min, min_Npoints )
    #
    return l1 + l2 + l3 + l4



colors_red  = ['#550000', '#994444', '#ee7777' ]
colors_blue = ['#2A2871', '#4742E7', '#D7D6FF' ]


def draw_likelihood_2D( fig, plot, vector_weights, matrix_parameter, ranges=[], grid_size=100, draw='likelihood contour', **kwargs ):
    if len(ranges)==0:
        ranges = np.zeros( (2,2) )
        ranges[:,0] = np.amin(matrix_parameter, axis=0)
        ranges[:,1] = np.amax(matrix_parameter, axis=0)
    #
    cmap                ='magma_r'
    smoothing           =False
    rel_sigma_smoothing = 0.025
    for key in kwargs:
        if key=='cmap':
            cmap = kwargs['cmap']
        if key=='smoothing':
            smoothing = kwargs['smoothing']
        if key=='rel_sigma_smoothing':
            rel_sigma_smoothing = kwargs['rel_sigma_smoothing']
    #
    dx = (ranges[0,1]-ranges[0,0])/grid_size
    dy = (ranges[1,1]-ranges[1,0])/grid_size
    #
    g_x = np.arange(  ranges[0,0], ranges[0,1]+dx/2., dx  )
    g_y = np.arange(  ranges[1,0], ranges[1,1]+dy/2., dy  )
    #
    m_x, m_y = np.meshgrid(g_x, g_y, indexing='ij')
    #
    L = np.zeros( (grid_size,grid_size) )
    for i, w in enumerate( vector_weights ):
        b_x = int(0.999999*(matrix_parameter[i,0]-ranges[0,0])/dx)
        b_y = int(0.999999*(matrix_parameter[i,1]-ranges[1,0])/dy)
        if b_x>=0 and b_y>=0 and b_x<grid_size and b_y<grid_size:
            L[b_x,b_y] += w
    #
    # smoothing (Gaussian)
    if smoothing:
        Ls = np.zeros( (grid_size,grid_size) )
        sigma   = rel_sigma_smoothing
        sigma_i = sigma*grid_size
        #
        ip, jp = np.meshgrid( np.arange(grid_size), np.arange(grid_size), indexing='ij'  )
        for i in range(grid_size):
            #print( i )
            for j in range(grid_size):
                Ls[i,j] = ( L[ip, jp] * np.exp( - ( (i-ip)**2 + (j-jp)**2 )/2./sigma_i**2 ) ).sum()
        Ls *= 1./(2*math.pi*sigma_i**2)
        L = Ls
    #
    #print(L.sum())
    L /= L.sum()
    #
    # get levels:
    L_sort_1D = np.sort( L.flatten() )
    levels = [0]
    sum = 0
    for i in range(len(L_sort_1D)):
        sum += L_sort_1D[i]
        if sum > 1. - 0.6827:
            levels.append(0.5*(L_sort_1D[i]+L_sort_1D[i-1]))
            break
        elif sum > 1. - 0.9545 and len(levels)==2:
            levels.append(0.5*(L_sort_1D[i]+L_sort_1D[i-1]))
        elif sum > 1. - 0.9973 and len(levels)==1:
            levels.append(0.5*(L_sort_1D[i]+L_sort_1D[i-1]))
    levels.append(1)
    #
    if not 'colors' in kwargs:
        kwargs['colors'] = colors_blue
    if not 'alpha' in kwargs:
        kwargs['alpha'] = 1.0
    if not 'ec'     in kwargs:
        kwargs['ec'] = 'black'
    if not 'lw' in kwargs:
        kwargs['lw'] = 1
    if not 'zorder' in kwargs:
        kwargs['zorder'] = 1
        
    if 'draw_sigmas' in kwargs:
        ds = kwargs['draw_sigmas']
        if len(ds)==2:
            if ds[0]==1 and ds[1]==2:
                levels = np.array(levels)
                #print(levels)
                kwargs['colors'] = kwargs['colors'][:2]
                levels = levels[ [0,2,3,4] ]
                #print(levels)
                
    #
    if 'contour'  in draw:
        #print(kwargs['ec'])
        plot.contour(m_x[:-1,:-1]+dx/2., m_y[:-1,:-1]+dy/2., L, levels[1:-1], colors=kwargs['ec'], linewidths=kwargs['lw'], zorder=10+kwargs['zorder']   )
        #cs = plot.contour(m_x[:-1,:-1]+dx/2., m_y[:-1,:-1]+dy/2., L, levels[1:-1], alpha=0. )
        #for col in cs.collections:
        #    p = col.get_paths()[0]
        #    v = p.vertices
        #    x = v[:,0]
        #    y = v[:,1]
        #    plot.plot( x, y, color=kwargs['ec'], lw=kwargs['lw'], zorder=10+kwargs['zorder']  )
    if 'contourf' in draw:
        #print(kwargs['alpha'])
        plot.contourf(m_x[:-1,:-1]+dx/2., m_y[:-1,:-1]+dy/2., L, levels[1:], colors=(kwargs['colors'])[::-1], alpha=kwargs['alpha'], zorder=kwargs['zorder'] )
    if 'likelihood' in draw:
        cmap = plt.get_cmap(cmap)
        norm = mpl.colors.Normalize(vmin=0, vmax=np.amax(L))
        plot.pcolormesh(  m_x, m_y, L, norm=norm, cmap=cmap, zorder=-10+kwargs['zorder']  )


def draw_likelihood_2D_adaptive_grid( fig, plot, vector_weights, matrix_parameter, **kwargs ):
    plot.set_xlim(  ( np.amin(matrix_parameter[:,0]),np.amax(matrix_parameter[:,0]) )  )
    plot.set_ylim(  ( np.amin(matrix_parameter[:,1]),np.amax(matrix_parameter[:,1]) )  )
    if not 'r_min' in kwargs:
        kwargs['r_min'] = 3
    if not 'r_max' in kwargs:
        kwargs['r_max'] = 6
    if not 'min_Npoints' in kwargs:
        kwargs['min_Npoints'] = 15
    bins = np.array( get_bayes_in_adaptive_grid(matrix_parameter[:,0], matrix_parameter[:,1], vector_weights, r_min=kwargs['r_min'],r_max=kwargs['r_max'],min_Npoints=kwargs['min_Npoints']) )
    cmap='magma_r'
    alpha  = 1.0
    for key in kwargs:
        if key=='cmap':
            cmap = kwargs[key]
        if key=='alpha':
            alpha  = kwargs[key]
    cmap = plt.get_cmap(cmap)
    draw_grid = False
    if 'draw_grid' in kwargs:
        draw_grid = kwargs['draw_grid']
    lw = 0.5
    if 'lw' in kwargs:
        lw = kwargs['lw']
    Lmax = np.amax(bins[:,4])
    for bin in bins:
        c = bin[4]/Lmax
        plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=cmap(c), alpha=alpha) )
        if draw_grid:
            plot.plot( (bin[0],bin[0],bin[1]), (bin[3],bin[2],bin[2]), color='black', lw=lw  )


def draw_contours_2D_1to3sigma_adaptive_grid( fig, plot, vector_weights, matrix_parameter, **kwargs ):
    plot.set_xlim(  ( np.amin(matrix_parameter[:,0]),np.amax(matrix_parameter[:,0]) )  )
    plot.set_ylim(  ( np.amin(matrix_parameter[:,1]),np.amax(matrix_parameter[:,1]) )  )
    if not 'r_min' in kwargs:
        kwargs['r_min'] = 3
    if not 'r_max' in kwargs:
        kwargs['r_max'] = 6
    if not 'min_Npoints' in kwargs:
        kwargs['min_Npoints'] = 15
    bins = np.array( get_bayes_in_adaptive_grid(matrix_parameter[:,0], matrix_parameter[:,1], vector_weights, r_min=kwargs['r_min'],r_max=kwargs['r_max'],min_Npoints=kwargs['min_Npoints']) )
    ind = np.argsort( bins[:,4] )
    bins[:,:] = bins[ind[::-1],:]
    colors = colors_red
    alpha  = 1.0
    for key in kwargs:
        if   key=='colors':
            colors = kwargs[key]
        elif key=='alpha':
            alpha  = kwargs[key]
    draw_grid = False
    if 'draw_grid' in kwargs:
        draw_grid = kwargs['draw_grid']
    lw = 0.5
    if 'lw' in kwargs:
        lw = kwargs['lw']
    Lmax = np.amax(bins[:,4])
    s = 0
    for bin in bins:
        s += bin[4]*bin[5]
        if s < 0.6827:
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=colors[0], alpha=alpha) )
        elif s < 0.9545:
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=colors[1], alpha=alpha) )
        elif s < 0.9973:
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=colors[2], alpha=alpha) )
        if draw_grid:
            plot.plot( (bin[0],bin[0],bin[1]), (bin[3],bin[2],bin[2]), color='black', lw=lw  )

def draw_histogram_1D(fig, plot, vector_weights, vector_parameter, **kwargs):
    if not 'n_histbins' in kwargs:
        kwargs['n_histbins'] = 50
    if not 'color' in kwargs:
        kwargs['color'] = 'black'
    if not 'norm_1' in kwargs:
        kwargs['norm_1'] = False
    if not 'smoothing' in kwargs:
        kwargs['smoothing'] = 0
    #
    p_from = np.amin(vector_parameter)
    p_to   = np.amax(vector_parameter)
    p_step = (p_to - p_from)/(kwargs['n_histbins'])
    #
    hist_bins  = np.arange(p_from,p_to+0.5*p_step,p_step)
    hist_means = 0.5 * ( hist_bins[1:] + hist_bins[:-1] )
    hist_values= np.zeros(len(hist_means))
    #
    for j, p in enumerate(vector_parameter):
        w = vector_weights[j]
        bin  = int(kwargs['n_histbins']*(p-p_from)/(p_to-p_from))
        if bin>=0 and bin<kwargs['n_histbins']:
            hist_values[bin] += w
    #
    smoothing = kwargs['smoothing']
    ## We multiply with a Gaussian kernal of width: <smoothing> * bin_width
    if smoothing > 0 :
        new_values     = np.copy(hist_values)
        new_values     = np.log(new_values+1e-90)
        new_values_II  = np.copy(new_values)
        for i,v in enumerate(new_values):
            new_values_II[i] = 0
            for j in range(-10, len(hist_values)+10):
                jj = j
                if j < 0:
                    jj = 0
                if j >= len(hist_values):
                    jj = len(hist_values) - 1
                new_values_II[i] += 1./np.sqrt( 2* math.pi * smoothing**2 ) * np.exp( - (i-j)**2/2./smoothing**2 ) * new_values[jj]
        hist_values = np.exp(new_values_II)
    #
    p_x = np.zeros(2*len(hist_means))
    p_y = np.copy(p_x)
    #
    p_x[0::2] = hist_bins[ :-1]
    p_x[1::2] = hist_bins[ 1: ]
    p_y[0::2] = hist_values[:]
    p_y[1::2] = hist_values[:]
    #
    if not 'style' in kwargs:
        kwargs['style'] = 'hist'
    norm=1.
    if kwargs['norm_1']:
        norm = np.amax(hist_values)
    if   kwargs['style'] == 'hist':
        plot.plot(p_x,        p_y/norm,         color=kwargs['color'], lw=1.5)
    elif kwargs['style'] == 'line':
        plot.plot(hist_means, hist_values/norm, color=kwargs['color'], lw=1.5)
    plot.set_xlim( (p_from,p_to) )
    plot.set_ylim( (0, 1.2*np.amax(hist_values)/norm ) )
