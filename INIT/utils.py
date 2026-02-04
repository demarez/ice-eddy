import numpy as np


def smooth_tanh(x,w=None,x_0=None,μ=1):
    '''
    x   [array] : Space to evaluate function
    w  [float] : Width of smooth transition
                  Default (None), i.e. 10% of the max value of xarray
    x_0 [float] : Center of smooth transition
                  Default (None), i.e. half of the max value of xarray
    '''
    if type(x_0) != np.ndarray or type(w) == np.ndarray:
        if x_0 == None:
            x_0=np.max(x)/2
        if np.array(w).any() == None:
            w = 0.1*np.max(x)
    else:
        x=np.ones((len(x_0),*x.shape)) * x
        x_0_tmp = x_0.copy()
        for ii in range(len(x.shape) - 1):
            x_0_tmp = np.expand_dims(x_0_tmp, axis=1)
        x_0 = x_0_tmp
    
    s = 1/2+1/2*np.tanh((x-x_0)/(μ/(np.pi*(2/w))))
    return s


def compute_depth_prof(ngrid_z, max_depth):
    '''
    Function to compute non regular depth grid
    '''
    dept_w=[0]
    dept_t=[0.5]
    
    nz = ngrid_z

    tah_max=(((ngrid_z) - 0.8*ngrid_z))/200
    scale = (ngrid_z-1)+(1/2+1/2*np.tanh(2*np.pi*tah_max))*(max_depth-(ngrid_z-1))/(ngrid_z-1)*(ngrid_z-1)

    for jk in range(1,nz):
        tah_val=(((jk+1) - 0.8*ngrid_z))/200

        dept=jk+(1/2+1/2*np.tanh(2*np.pi*tah_val))*(max_depth-(ngrid_z-1))/(ngrid_z-1)*jk
        z_tanh_z_norm = (dept/scale)*max_depth
        dept_w.append(z_tanh_z_norm)
        dept_t.append(z_tanh_z_norm + (z_tanh_z_norm-dept_w[-2])/2)
        
    return dept_w,dept_t
