import  numpy                   as      np
import  math

import  matplotlib              as      mpl
mpl.use('Agg')
import  matplotlib.pyplot       as      plt
import  matplotlib.patches      as      mpatches
from matplotlib.patches         import  Rectangle


colors_red  = ['#880000', '#bb4444', '#ff7777' ]
colors_blue = ['#2A2871', '#4742E7', '#D7D6FF' ]
colors_green= ['#005500', '#009900', '#00ff00' ]
colors_amber= ['#ED872D', '#FF7E00', '#FFA700' ]
colors_grey = ['#444444', '#888888', '#aaaaaa' ]
colors_cyan = ['#24B3B8', '#6FE4E8', '#9DF0F3' ]
colors_mag  = ['#71286D', '#BA10B1', '#DC96D8' ]


def get_chi2_profile_in_adaptive_grid__old( vec_x, vec_y, vec_chi2, range_x=None, range_y=None, r=0, r_max=6, r_min=3, min_Npoints=15):
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
            chi2=1e90
            if len(vec_chi2):
                chi2 = np.amin(vec_chi2)
            return [[ x1, x3, y1, y3, chi2]]
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
        chi2 = vec_chi2[i]
        if x < x1 or x > x3 or y < y1 or y>y3:
            continue
        if x < x2:
            if y < y2:
                ll_vx.append(x)
                ll_vy.append(y)
                ll_chi2.append(chi2)
            else:
                ul_vx.append(x)
                ul_vy.append(y)
                ul_chi2.append(chi2)
        else:
            if y < y2:
                lr_vx.append(x)
                lr_vy.append(y)
                lr_chi2.append(chi2)
            else:
                ur_vx.append(x)
                ur_vy.append(y)
                ur_chi2.append(chi2)
    #
    l1 = get_chi2_profile_in_adaptive_grid( ll_vx, ll_vy, ll_chi2, [x1,x2], [y1,y2], r+1, r_max, r_min, min_Npoints )
    l2 = get_chi2_profile_in_adaptive_grid( ul_vx, ul_vy, ul_chi2, [x1,x2], [y2,y3], r+1, r_max, r_min, min_Npoints )
    l3 = get_chi2_profile_in_adaptive_grid( lr_vx, lr_vy, lr_chi2, [x2,x3], [y1,y2], r+1, r_max, r_min, min_Npoints )
    l4 = get_chi2_profile_in_adaptive_grid( ur_vx, ur_vy, ur_chi2, [x2,x3], [y2,y3], r+1, r_max, r_min, min_Npoints )
    #
    return l1 + l2 + l3 + l4


def get_chi2_profile_in_adaptive_grid(vec_x, vec_y, vec_chi2, range_x=None, range_y=None, A_min=1e-3, A_max=1e-2, max_ratio=4, min_Npoints=15):
    
    # determine range
    if range_x==None:
        range_x = [np.amin(vec_x), np.amax(vec_x) ]
    if range_y==None:
        range_y = [np.amin(vec_y), np.amax(vec_y) ]
    x_l = range_x[0]
    x_u = range_x[1]
    y_l = range_y[0]
    y_u = range_y[1]

    my_vec_x    = []
    my_vec_y    = []
    my_vec_chi2 = []

    for i in range(len(vec_x)):
        x = vec_x   [i]
        y = vec_y   [i]
        c = vec_chi2[i]
        if x<=x_u and x>=x_l and y<=y_u and y>=y_l:
            my_vec_x   .append(x)
            my_vec_y   .append(y)
            my_vec_chi2.append(c)

    my_vec_x     = np.array( my_vec_x    )
    my_vec_y     = np.array( my_vec_y    )
    my_vec_chi2  = np.array( my_vec_chi2 )

    t_vec_x = (my_vec_x-x_l)/(x_u-x_l)
    t_vec_y = (my_vec_y-y_l)/(y_u-y_l)
    res = _get_chi2_profile_in_adaptive_grid__parameteter_between_0_1(t_vec_x, t_vec_y, my_vec_chi2, range_x=[0.,1.], range_y=[0.,1.], A_min=A_min, A_max=A_max, max_ratio=max_ratio, min_Npoints=min_Npoints)

    res = np.array(res)

    res[:,0] = x_l+(x_u-x_l)*res[:,0]
    res[:,1] = x_l+(x_u-x_l)*res[:,1]
    res[:,2] = y_l+(y_u-y_l)*res[:,2]
    res[:,3] = y_l+(y_u-y_l)*res[:,3]

    return res


def _get_chi2_profile_in_adaptive_grid__parameteter_between_0_1(vec_x, vec_y, vec_chi2, range_x=[0.,1.], range_y=[0.,1.], A_min=1e-3, A_max=1e-2, max_ratio=4, min_Npoints=50):

    x_l = range_x[0]
    x_u = range_x[1]
    y_l = range_y[0]
    y_u = range_y[1]
    A   = (x_u-x_l)*(y_u-y_l)
    
    if len(vec_x)==0:
        return [[ x_l, x_u, y_l, y_u, 1e90, A]]
    
    if A<=A_max:
        if (A<=A_min*max_ratio or len(vec_x)<min_Npoints):
            chi2 = np.amin(vec_chi2)
            return [[ x_l, x_u, y_l, y_u, chi2, A]]
    
    splitX = True
    if (y_u-y_l)>(x_u-x_l):
        splitX = False

    if splitX:
        ind = np.argsort( vec_x )
    else:
        ind = np.argsort( vec_y )
    vec_x   [:] = vec_x   [ind]
    vec_y   [:] = vec_y   [ind]
    vec_chi2[:] = vec_chi2[ind]

    i_cut = int( (len(vec_x)-1)/2 )
    if splitX:
        x_cut_min = x_l + (  1./max_ratio) * (x_u-x_l)
        x_cut_max = x_l + (1-1./max_ratio) * (x_u-x_l)
        x_cut = vec_x[i_cut]
        if x_cut<x_cut_min:
            x_cut = x_cut_min
        if x_cut>x_cut_max:
            x_cut = x_cut_max
        i_cut = 0
        for i in range(len(vec_x)):
            if vec_x[i]<x_cut:
                i_cut+=1
        l1 = _get_chi2_profile_in_adaptive_grid__parameteter_between_0_1( vec_x[:i_cut], vec_y[:i_cut], vec_chi2[:i_cut], [x_l,x_cut], [y_l,y_u], A_min, A_max, max_ratio, min_Npoints )
        l2 = _get_chi2_profile_in_adaptive_grid__parameteter_between_0_1( vec_x[i_cut:], vec_y[i_cut:], vec_chi2[i_cut:], [x_cut,x_u], [y_l,y_u], A_min, A_max, max_ratio, min_Npoints )
        return l1 + l2
    else:
        y_cut_min = y_l + (  1./max_ratio) * (y_u-y_l)
        y_cut_max = y_l + (1-1./max_ratio) * (y_u-y_l)
        y_cut = vec_y[i_cut]
        if y_cut<y_cut_min:
            y_cut = y_cut_min
        if y_cut>y_cut_max:
            y_cut = y_cut_max
        i_cut = 0
        for i in range(len(vec_y)):
            if vec_y[i]<y_cut:
                i_cut+=1
        l1 = _get_chi2_profile_in_adaptive_grid__parameteter_between_0_1( vec_x[:i_cut], vec_y[:i_cut], vec_chi2[:i_cut], [x_l,x_u], [y_l,y_cut], A_min, A_max, max_ratio, min_Npoints )
        l2 = _get_chi2_profile_in_adaptive_grid__parameteter_between_0_1( vec_x[i_cut:], vec_y[i_cut:], vec_chi2[i_cut:], [x_l,x_u], [y_cut,y_u], A_min, A_max, max_ratio, min_Npoints )
        return l1 + l2


#### Fuctions to draw a contour arround points

def _xy_to_phi(x,y):
    d=0
    if y < 0:
        d=2*math.pi
    return np.arctan2( y,x )+d
xy_to_phi = np.vectorize(_xy_to_phi)
def _to_positive_angle(phi):
    if phi < 0:
        return _to_positive_angle(phi + 2*math.pi)
    return phi
to_positive_angle = np.vectorize(_to_positive_angle)
def get_next_point(x_vec, y_vec, index, angle, relative_interpolation_lenght=0.3 ):
    mp1_vec =  np.cos(angle) * x_vec  + np.sin(angle) * y_vec
    mp2_vec = -np.sin(angle) * x_vec  + np.cos(angle) * y_vec
    angles = xy_to_phi( mp1_vec-mp1_vec[index],mp2_vec-mp2_vec[index] )
    diffs  = np.sqrt(  (mp1_vec-mp1_vec[index])**2 + (mp2_vec-mp2_vec[index])**2 )
    for i, d in enumerate(diffs):
        if d>relative_interpolation_lenght                           or d<0.001:
        #if d>relative_interpolation_lenght/(0.001+np.sin(angles[i])) or d<0.001:
            angles[i] = 100
    new_index = np.argmin(angles)
    new_angle = xy_to_phi( x_vec[new_index]-x_vec[index],y_vec[new_index]-y_vec[index] )
    return new_index, to_positive_angle(new_angle - math.pi*4./8)
def get_contour( x_vec, y_vec,  relative_interpolation_lenght=0.3 ):
    p1_vec = ( x_vec - np.amin(x_vec) )/(np.amax(x_vec)-np.amin(x_vec))
    p2_vec = ( y_vec - np.amin(y_vec) )/(np.amax(y_vec)-np.amin(y_vec))
    boarder_points = []
    index_start = np.argmax(p1_vec)
    boarder_points.append(index_start)
    current_index = index_start
    current_angel = 0
    n = 100
    for i in range(n):
        i,current_angel = get_next_point(p1_vec, p2_vec, boarder_points[-1], current_angel, relative_interpolation_lenght)
        #print(i)
        if i in boarder_points:
            boarder_points.append(i)
            break
        boarder_points.append(i)
    return x_vec[boarder_points], y_vec[boarder_points], np.array(boarder_points)
import sys
sys.setrecursionlimit(3000)
def recursive_cluster_finder( x_vec, y_vec, cluster_vec=[], distSq=0.1**2, i=0, c=1  ):      # put this to cpp
    if len(cluster_vec)==0:
        x_vec=np.copy(x_vec)
        y_vec=np.copy(y_vec)
        x_vec = ( x_vec - np.amin(x_vec) )/(np.amax(x_vec)-np.amin(x_vec))
        y_vec = ( y_vec - np.amin(y_vec) )/(np.amax(y_vec)-np.amin(y_vec))
        
        cluster_vec = np.ones(len(x_vec),dtype=int)*1001
        cluster_vec[i] = c
    for j in range(len(x_vec)):
        if j==i:
            continue
        if (x_vec[i]-x_vec[j])**2 + (y_vec[i]-y_vec[j])**2 < distSq:
            if not cluster_vec[j]==c:
                cluster_vec[j]=c
                recursive_cluster_finder( x_vec, y_vec, cluster_vec, distSq, i=j, c=c  )
    if i==0:
        for j in range(len(x_vec)):
            if cluster_vec[j]>1000:
                c=c+1
                recursive_cluster_finder( x_vec, y_vec, cluster_vec, distSq, i=j, c=c  )
        return cluster_vec
    return

from matplotlib.path import Path
import matplotlib.patches as patches

def xy_to_path(x,y,):
    verts = []
    codes = [Path.MOVETO]
    for i,_x in enumerate(x):
        _y = y[i]
        verts.append((_x,_y))
        if i>0:
            codes.append(Path.LINETO)
    verts.append(verts[0])
    codes.append(Path.CLOSEPOLY)

    path = Path(verts, codes)
    return path

# one might consider to apply a cluster finding first: recursive_cluster_finder
def draw_cluster_contour( fig, plot, vector_chi2, matrix_parameter, **kwargs ):
    chi2_min = vector_chi2[0]
    for i, chi2 in enumerate(vector_chi2):
        if   chi2-chi2_min < 2.3:
            ind_1s_2D = i
        elif chi2-chi2_min < 6.18:
            ind_2s_2D = i
        elif chi2-chi2_min < 11.83:
            ind_3s_2D = i
        else:
            break
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

    if not 'relative_interpolation_lenght' in kwargs:
        kwargs['relative_interpolation_lenght'] = 0.3

    # implemennt cluster finder (optional) 

    x_1s, y_1s, _ = get_contour(matrix_parameter[:ind_1s_2D, 0], matrix_parameter[:ind_1s_2D, 1],relative_interpolation_lenght=kwargs['relative_interpolation_lenght'])
    x_2s, y_2s, _ = get_contour(matrix_parameter[:ind_2s_2D, 0], matrix_parameter[:ind_2s_2D, 1],relative_interpolation_lenght=kwargs['relative_interpolation_lenght'])
    x_3s, y_3s, _ = get_contour(matrix_parameter[:ind_3s_2D, 0], matrix_parameter[:ind_3s_2D, 1],relative_interpolation_lenght=kwargs['relative_interpolation_lenght'])

    patch1 = patches.PathPatch(xy_to_path(x_1s[:-1], y_1s[:-1]), facecolor=kwargs['colors'][0], lw=0, alpha=kwargs['alpha'], zorder=kwargs['zorder'] )
    patch2 = patches.PathPatch(xy_to_path(x_2s[:-1], y_2s[:-1]), facecolor=kwargs['colors'][1], lw=0, alpha=kwargs['alpha'], zorder=kwargs['zorder']  )
    patch3 = patches.PathPatch(xy_to_path(x_3s[:-1], y_3s[:-1]), facecolor=kwargs['colors'][2], lw=0, alpha=kwargs['alpha'], zorder=kwargs['zorder'] )

    plot.add_patch(patch3)
    plot.add_patch(patch2)
    plot.add_patch(patch1)

    plot.plot( x_1s, y_1s, color=kwargs['ec'], lw=kwargs['lw'], zorder=(kwargs['zorder']+100) )
    plot.plot( x_2s, y_2s, color=kwargs['ec'], lw=kwargs['lw'], zorder=(kwargs['zorder']+100) )
    plot.plot( x_3s, y_3s, color=kwargs['ec'], lw=kwargs['lw'], zorder=(kwargs['zorder']+100) )


def draw_best_fitpoint( fig, plot, vector_chi2, matrix_parameter, **kwargs ):
    if not 'alpha' in kwargs:
        kwargs['alpha'] = 1.0
    if not 'ec'     in kwargs:
        kwargs['ec'] = 'black'
    if not 'fc'     in kwargs:
        kwargs['fc'] = 'white'
    if not 'lw' in kwargs:
        kwargs['lw'] = 1
    if not 'markersize' in kwargs:
        kwargs['markersize'] = 8
    if not 'zorder' in kwargs:
        kwargs['zorder'] = 1000
    if not 'fmt' in kwargs:
        kwargs['fmt'] = 'o'
    if not 'vector_positions' in kwargs:
        kwargs['vector_positions'] = [0,1]
    xerr = 0
    yerr = 0
    if 'err' in kwargs:
        xerr = (kwargs['err'])[kwargs['vector_positions'][0]]
        yerr = (kwargs['err'])[kwargs['vector_positions'][1]]

    x = matrix_parameter[0,0]
    y = matrix_parameter[0,1]

    plot.errorbar( x, y, xerr=[xerr], yerr=[yerr], fmt=kwargs['fmt'], color=kwargs['ec'], mfc=kwargs['fc'], lw=kwargs['lw'], markersize=kwargs['markersize'], zorder=kwargs['zorder'] )


def draw_scatter( fig, plot, vector_chi2, matrix_parameter, **kwargs ):
    lw=0
    range_delta_chi2 = 20
    cmap='magma'
    s=5
    kwargs_forward = {}
    for key in kwargs:
        if key=='s':
            s = kwargs['s']
        elif key=='lw':
            lw = kwargs['lw']
        elif key=='range_delta_chi2':
            range_delta_chi2 = kwargs['range_delta_chi2']
        elif key=='cmap':
            cmap = kwargs['cmap']
        elif key=='vector_positions':
            continue
        else:
            kwargs_forward[key] = kwargs[key]
    cmap = mpl.cm.get_cmap(cmap)
    im = plot.scatter( matrix_parameter[::-1,0],matrix_parameter[::-1,1], c=vector_chi2[::-1], lw=lw, vmin=vector_chi2[0], vmax=vector_chi2[0]+range_delta_chi2, cmap=cmap, s=s, **kwargs_forward )
    return im 

def draw_scatter_1to3sigma( fig, plot, vector_chi2, matrix_parameter, **kwargs ):
    lw=0
    s=5
    colors = colors_red
    alpha  = 1.0
    zorder = 1
    kwargs_forward = {}
    for key in kwargs:
        if   key=='s':
            s = kwargs[key]
        elif key=='lw':
            lw = kwargs[key]
        elif key=='colors':
            colors = kwargs[key]
        elif key=='alpha':
            alpha  = kwargs[key]
        elif key=='zorder':
            zorder = kwargs[key]
        elif key=='vector_positions':
            continue
        else:
            kwargs_forward[key] = kwargs[key]
    #
    chi2_min = np.amin(vector_chi2)
    ## Find 1 to 3 sigma intervals
    ind_1s_2D = 0
    ind_2s_2D = 0
    ind_3s_2D = 0
    #
    for i, chi2 in enumerate(vector_chi2):
        if   chi2-chi2_min < 2.3:
            ind_1s_2D = i
        elif chi2-chi2_min < 6.18:
            ind_2s_2D = i
        elif chi2-chi2_min < 11.83:
            ind_3s_2D = i
        else:
            break
    plot.scatter( matrix_parameter[ind_2s_2D:ind_3s_2D,0],matrix_parameter[ind_2s_2D:ind_3s_2D,1], c=colors[2], lw=lw, s=s, zorder=zorder, **kwargs_forward )
    plot.scatter( matrix_parameter[ind_1s_2D:ind_2s_2D,0],matrix_parameter[ind_1s_2D:ind_2s_2D,1], c=colors[1], lw=lw, s=s, zorder=zorder, **kwargs_forward )
    plot.scatter( matrix_parameter[         :ind_1s_2D,0],matrix_parameter[         :ind_1s_2D,1], c=colors[0], lw=lw, s=s, zorder=zorder, **kwargs_forward )


def draw_contour_2D_1to3sigma_adaptive_grid( fig, plot, vector_chi2, matrix_parameter, **kwargs ):
    plot.set_xlim(  ( np.amin(matrix_parameter[:,0]),np.amax(matrix_parameter[:,0]) )  )
    plot.set_ylim(  ( np.amin(matrix_parameter[:,1]),np.amax(matrix_parameter[:,1]) )  )
    if not 'A_min' in kwargs:
        kwargs['A_min'] = 1e-3
    if not 'A_max' in kwargs:
        kwargs['A_max'] = 1e-2
    if not 'max_ratio' in kwargs:
        kwargs['max_ratio'] = 4
    if not 'min_Npoints' in kwargs:
        kwargs['min_Npoints'] = 15
    bins = np.array( get_chi2_profile_in_adaptive_grid(matrix_parameter[:,0], matrix_parameter[:,1], vector_chi2, A_min=kwargs['A_min'],A_max=kwargs['A_max'],max_ratio=kwargs['max_ratio'],min_Npoints=kwargs['min_Npoints']) )
    chi2_min = np.amin(vector_chi2)
    colors = colors_red
    alpha  = 1.0
    for key in kwargs:
        if key=='colors':
            colors = kwargs[key]
        elif key=='alpha':
            alpha  = kwargs[key]
    draw_grid = False
    if 'draw_grid' in kwargs:
        draw_grid = kwargs['draw_grid']
    lw = 0.5
    if 'lw' in kwargs:
        lw = kwargs['lw']
    for bin in bins:
        bin_min_chi2 = bin[4]-chi2_min
        if bin_min_chi2 <= 2.3:
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=colors[0], alpha=alpha) )
        if bin_min_chi2<=6.18 and bin_min_chi2>2.3:
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=colors[1], alpha=alpha) )
        if bin_min_chi2<=11.83 and bin_min_chi2>6.18 :
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=colors[2], alpha=alpha) )
        if draw_grid:
            plot.plot( (bin[0],bin[0],bin[1]), (bin[3],bin[2],bin[2]), color='black', lw=lw  )

def draw_chi2_profile_adaptive_grid( fig, plot, vector_chi2, matrix_parameter, **kwargs ):
    
    range_delta_chi2 = 20
    cmap='magma'
    for key in kwargs:
        if key=='cmap':
            cmap = kwargs['cmap']
        elif key=='range_delta_chi2':
            range_delta_chi2 = kwargs['range_delta_chi2']
    cmap = mpl.cm.get_cmap('magma')
    
    plot.set_xlim(  ( np.amin(matrix_parameter[:,0]),np.amax(matrix_parameter[:,0]) )  )
    plot.set_ylim(  ( np.amin(matrix_parameter[:,1]),np.amax(matrix_parameter[:,1]) )  )
    if not 'A_min' in kwargs:
        kwargs['A_min'] = 1e-3
    if not 'A_max' in kwargs:
        kwargs['A_max'] = 1e-2
    if not 'max_ratio' in kwargs:
        kwargs['max_ratio'] = 4
    if not 'min_Npoints' in kwargs:
        kwargs['min_Npoints'] = 15
    bins = np.array( get_chi2_profile_in_adaptive_grid(matrix_parameter[:,0], matrix_parameter[:,1], vector_chi2, A_min=kwargs['A_min'],A_max=kwargs['A_max'],max_ratio=kwargs['max_ratio'],min_Npoints=kwargs['min_Npoints']) )
    chi2_min = np.amin(vector_chi2)
    alpha  = 1.0
    for key in kwargs:
        if key=='alpha':
            alpha  = kwargs[key]
    draw_grid = False
    if 'draw_grid' in kwargs:
        draw_grid = kwargs['draw_grid']
    lw = 0.5
    if 'lw' in kwargs:
        lw = kwargs['lw']
    for bin in bins:
        bin_min_chi2 = bin[4]-chi2_min
        if bin_min_chi2 <= range_delta_chi2:
            c = bin_min_chi2/range_delta_chi2
            plot.add_patch( Rectangle( (bin[0],bin[2]), bin[1]-bin[0], bin[3]-bin[2], color=cmap(c), alpha=alpha) )
        if draw_grid:
            plot.plot( (bin[0],bin[0],bin[1]), (bin[3],bin[2],bin[2]), color='black', lw=lw  )
    im = plot.scatter( matrix_parameter[::-1,0], matrix_parameter[::-1,1], c=vector_chi2[::-1],  lw=0, vmin=chi2_min, vmax=chi2_min+range_delta_chi2, s=0, cmap=cmap )
    return im
    #draw_split(rx, ry, vx, vy, vchi2)
    #plt.subplots_adjust(left=0.2, right=0.7, top=0.9, bottom=0.15)
    #cbar_ax = fig.add_axes([0.75,0.15,0.05,0.75])
    #o = fig.colorbar( mappable=im, cax=cbar_ax  )
    #o.set_label(r'$\chi^{2}}$', rotation=0, labelpad=5, fontsize=20)
    #plot.set_xlim( rx )
    #plot.set_ylim( ry )
    #plt.savefig(wdir+'/results/'+run_name+'/chi2_adaptiveGrid__%s_%s.png'  %(list_fit_par[i], list_fit_par[j]))



def draw_histogram_1D(fig, plot, vector_chi2, vector_parameter, **kwargs):
    if not 'n_histbins' in kwargs:
        kwargs['n_histbins'] = 30
    if not 'color' in kwargs:
        kwargs['color'] = 'black'
    if not 'range_delta_chi2' in kwargs:
        kwargs['range_delta_chi2'] = 10
    if not 'plot_Delta_logLik' in kwargs:
        kwargs['plot_Delta_logLik'] = False
    if not 'smoothing' in kwargs:
        kwargs['smoothing'] = 0
    #
    p_from = np.amin(vector_parameter)
    p_to   = np.amax(vector_parameter)
    p_step = (p_to - p_from)/(kwargs['n_histbins'])
    if p_step==0.:
        return
    #
    chi2_min = np.amin(vector_chi2)
    hist_bins  = np.arange(p_from,p_to+0.5*p_step,p_step)
    hist_means = 0.5 * ( hist_bins[1:] + hist_bins[:-1] )
    hist_values= np.ones(len(hist_means))*(kwargs['range_delta_chi2']+10+chi2_min)
    #
    for j, p in enumerate(vector_parameter):
        chi2 = vector_chi2[j]
        bin  = int(kwargs['n_histbins']*(p-p_from)/(p_to-p_from))
        if bin>=0 and bin<kwargs['n_histbins']:
            hist_values[bin] = np.minimum( hist_values[bin], chi2 )
    #
    smoothing = kwargs['smoothing']
    ## We multiply with a Gaussian kernal with a width of <smoothing> * bin_width
    if smoothing > 0 :
        new_values     = np.copy(hist_values)
        new_values     = new_values - np.amin(hist_values) + 10.
        new_values     = np.log(new_values)
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
        hist_values = np.exp(new_values_II) - np.amin(np.exp(new_values_II)) + np.amin(hist_values)
    #
    p_x = np.zeros(2*len(hist_means))
    p_y = np.copy(p_x)
    p_x[0::2] = hist_bins[ :-1]
    p_x[1::2] = hist_bins[ 1: ]
    p_y[0::2] = hist_values[:]
    p_y[1::2] = hist_values[:]
    #
    if not 'style' in kwargs:
        kwargs['style'] = 'hist'
    d = 0.
    if kwargs['plot_Delta_logLik']:
        d = -chi2_min
    if   kwargs['style'] == 'hist':
        plot.plot(p_x,        p_y+d,         color=kwargs['color'], lw=1.5)
    elif kwargs['style'] == 'line':
        plot.plot(hist_means, hist_values+d, color=kwargs['color'], lw=1.5)
    #plot.set_xlim( (p_from,p_to) )
    plot.set_ylim( (chi2_min+d,chi2_min+d+kwargs['range_delta_chi2']) )


def get_bestfit_and_uncertainty(vector_chi2, vector_parameter, sigma=1):
    # Get min chiSq
    sorted     = np.argsort(vector_chi2)

    chi2_best  = vector_chi2      [sorted[0]]
    p_best     = vector_parameter [sorted[0]]
    
    p_u   = p_best
    p_l   = p_best
    chi_u = chi2_best
    chi_l = chi2_best
    for i in sorted:
        #print(vector_chi2[i] - chi2_best)
        if vector_chi2[i] - chi2_best > sigma**2:
            break
        if p_l > vector_parameter[i]:
            p_l   = vector_parameter[i]
            chi_l = vector_chi2[i]
        if p_u < vector_parameter[i]:
            p_u = vector_parameter[i]
            chi_u = vector_chi2[i]

    p_le   = p_l
    p_ue   = p_u
    for i in sorted:
        if vector_chi2[i] - chi2_best > (sigma+1)**2:
            break
        if p_l > vector_parameter[i]:
            m_l  = ( chi_l - vector_chi2[i] )/(p_l - vector_parameter[i])
            b_l  = ( chi_l * vector_parameter[i] - vector_chi2[i] * p_l )/(vector_parameter[i] - p_l)
            p_ll = ( chi2_best + sigma**2 - b_l )/m_l
            p_le = np.minimum(p_ll, p_le)
        if p_u < vector_parameter[i]:
            m_u  = ( chi_u - vector_chi2[i] )/(p_u - vector_parameter[i])
            b_u  = ( chi_u * vector_parameter[i] - vector_chi2[i] * p_u )/(vector_parameter[i] - p_u)
            p_uu = ( chi2_best + sigma**2 - b_u )/m_u
            p_ue = np.maximum(p_uu, p_ue)

    p_ue = np.minimum(np.amax(vector_parameter), p_ue)
    p_le = np.maximum(np.amin(vector_parameter), p_le)
#    print( '%8.2e  +%8.2e   -%8.2e'      % ( p_best, p_u-p_best,  p_best-p_l  )  )
#    print( '%8.2e  +%8.2e   -%8.2e  <--' % ( p_best, p_ue-p_best, p_best-p_le )  )
#    sys.exit(0)
    return p_best, p_ue-p_best, p_best-p_le
