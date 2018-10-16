
import numpy as np
cimport numpy as np
cimport cython

def particle_clump_mask_go(np.ndarray[np.float64_t, ndim=1] xpos ,
                           np.ndarray[np.float64_t, ndim=1] ypos ,
                           np.ndarray[np.float64_t, ndim=1] zpos,
                           np.ndarray[np.float64_t, ndim=1] grid_left,
                           np.ndarray[np.float64_t, ndim=1] grid_dx,
                           np.ndarray[np.int32_t, ndim=3] grid_mask
                          ):
    cdef int nparticles, nx, ny, nz
    nparticles = xpos.shape[0]
    cdef np.int64_t i,j,k,n
    cdef np.ndarray[np.int32_t, ndim=1] mask = np.zeros(nparticles, dtype='int32')

    for n in range( nparticles ):
        i = <int> ( (xpos[n] - grid_left[0])/grid_dx[0] )
        j = <int> ( (ypos[n] - grid_left[1])/grid_dx[1] )
        k = <int> ( (zpos[n] - grid_left[2])/grid_dx[2] )
        if grid_mask[i][j][k] == 1:
            mask[n] = 1
    return mask










#@particle_mask = particle_clump_mask( xpos, ypos, zpos, grid_left, grid_dx, clump_cut_mask)
