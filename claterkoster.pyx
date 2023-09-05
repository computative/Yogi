import cython
import numpy as np
cimport numpy as cnp
import utils
from libc.math cimport exp, pow
cnp.import_array()

def HGamma_m( int k, cnp.ndarray[cnp.float64_t, ndim=1] G, int N, double[:] param, cnp.ndarray[cnp.float64_t, ndim=2] R ):
    cdef double _ro = 2.35
    cdef double _rc = 3.8661
    cdef double _n = 1.9771
    cdef double _nc = 6.8702
    cdef double g = np.linalg.norm(G)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] r = np.zeros(3, dtype=np.float64)
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] HGamma = np.zeros( (N,N) , dtype=np.complex128 )
    # forst loop over offsite 
    cdef int i, j, L
    cdef double l, m, n, bowler
    cdef double IrI, x, y
    x = exp( -_n*pow(_rc,-_nc) )
    y = pow(_ro,_n)*exp(_n* pow( _ro/_rc, _nc) )
    L = len(R)
    for i in range( L ):
        for j in range(L):
            if i == j:
                continue
            r = R[j] + G - R[i]
            IrI = np.linalg.norm(r)
            e = r/IrI
            l,m,n = e[0], e[1], e[2]
            bowler = y*pow(x,(pow(IrI,_nc)))*pow(IrI,-_n)
            # what is bowler when r=g=0?
            HGamma[4*i:4*(i+1), 4*j:4*(j+1)] = np.array([
                [ bowler*param[2],      bowler*param[3]*l,                              bowler*param[3]*m,                              bowler*param[3]*n],
                [-bowler*param[3]*l,    bowler*param[4]*l*l + bowler*param[5]*(1-l*l ), bowler*m*l*(param[4]-param[5]),                 bowler*n*l*(param[4] - param[5])],
                [-bowler*param[3]*m,    bowler*l*m*(param[4] - param[5]),               bowler*param[4]*m*m + bowler*param[5]*(1-m*m),  bowler*n*m*(param[4] - param[5])],
                [-bowler*param[3]*n,    bowler*l*n*(param[4] - param[5]),               bowler*m*n*(param[4]-param[5]),                 bowler*param[4]*n*n + bowler*param[5]*(1-n*n)]])
    if g > utils.eps:
        for i in range(L):
            l, m, n = G/g
            # what is bowler when r=g=0?
            bowler = y*pow(x,(pow(g,_nc)))*pow(g,-_n)
            HGamma[4*i:4*(i+1),4*i:4*(i+1)] = np.array([
                [ bowler*param[2],      bowler*param[3]*l,                              bowler*param[3]*m,                              bowler*param[3]*n],
                [-bowler*param[3]*l,    bowler*param[4]*l*l + bowler*param[5]*(1-l*l),  bowler*m*l*(param[4]-param[5]),                 bowler*n*l*(param[4] - param[5])],
                [-bowler*param[3]*m,    bowler*l*m*(param[4] - param[5]),               bowler*param[4]*m*m + bowler*param[5]*(1-m*m),  bowler*n*m*(param[4] - param[5])],
                [-bowler*param[3]*n,    bowler*l*n*(param[4] - param[5]),               bowler*m*n*(param[4]-param[5]),                 bowler*param[4]*n*n + bowler*param[5]*(1-n*n)]])
    else:
        for i in range(L):
            HGamma[4*i:4*(i+1),4*i:4*(i+1)] = np.diag([param[0], param[1], param[1], param[1]]) 
    return HGamma

