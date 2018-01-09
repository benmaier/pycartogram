import numpy as np
import matplotlib.pyplot as pl

def fill_matrix(A,i,j,new_val=1.):
    index_list = [(i,j)]
    val = A[i,j]
    x_, y_ = A.shape

    while len(index_list) > 0:
        i,j = index_list.pop()
        A[i,j] = new_val
        if i+1 < x_ and A[i+1,j] == val:
            index_list.append((i+1,j))
        if j+1 < y_ and A[i,j+1] == val:
            index_list.append((i,j+1))
        if j-1 >= 0 and A[i,j-1] == val:
            index_list.append((i,j-1))
        if i-1 >= 0 and A[i-1,j] == val:
            index_list.append((i-1,j))            

    
        
def savefig_marginless(fn,fig,ax,**kwargs):
    ax.set_axis_off()
    fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    ax.margins(0,0)
    ax.xaxis.set_major_locator(pl.NullLocator())
    ax.yaxis.set_major_locator(pl.NullLocator())
    fig.savefig(fn, bbox_inches = 'tight',
    pad_inches = 0,**kwargs)

def scale(hist,
                     new_min=0.,
                     new_max=.9,
                     take_inverse=False,
                     replace_nan = True,
                     get_nan_min_max = False
                    ):
    if take_inverse:
        factor = -1
    else:
        factor = 1
    log_hist = factor * np.array(hist)
    nan_min = np.nanmin(log_hist)
    nan_max = np.nanmax(log_hist)
    if replace_nan:
        log_hist[np.isnan(log_hist)] = nan_min - 1
    min_ = nan_min - 1
    max_ = nan_max
    new_min = 0.
    new_max = 0.9
    intensity = lambda x: (x-min_)/ (max_-min_) * (new_max-new_min) + new_min
    if not get_nan_min_max:
        return log_hist, intensity
    else:
        return log_hist, intensity, nan_min, nan_max

def logify_and_scale(hist,
                     new_min=0.,
                     new_max=.9,
                     take_inverse=False,
                     replace_nan = True,
                     get_nan_min_max = False
                    ):
    if take_inverse:
        factor = -1
    else:
        factor = 1
    hist = np.array(hist)
    log_hist = np.array(hist)
    log_hist[hist>0.] = factor * np.log(hist[hist>0])
    log_hist[log_hist==0] = np.nan
    nan_min = np.nanmin(log_hist)
    nan_max = np.nanmax(log_hist)
    if replace_nan:
        log_hist[np.isnan(log_hist)] = nan_min - 1
    min_ = nan_min - 1
    max_ = nan_max
    new_min = 0.
    new_max = 0.9
    intensity = lambda x: (x-min_)/ (max_-min_) * (new_max-new_min) + new_min
    if not get_nan_min_max:
        return log_hist, intensity
    else:
        return log_hist, intensity, nan_min, nan_max

def is_iter(obj):
    try:
        some_object_iterator = iter(obj)
        return True
    except TypeError, te:
        return False

