import cython
import numpy as np
cimport numpy as cnp
import utils
from libc.math cimport exp, pow
cnp.import_array()

def HGamma( int k, cnp.ndarray[cnp.float64_t, ndim=1] G, int N, double[:] param, double[:] gsp_param, cnp.ndarray[cnp.float64_t, ndim=2] R ):
    cdef double _ro = gsp_param[0]
    cdef double _rc = gsp_param[1]
    cdef double _n = gsp_param[2]
    cdef double _nc = gsp_param[3]
    cdef double g = np.linalg.norm(G)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] r = np.zeros(3, dtype=np.float64)
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] HGamma = np.zeros( (N,N) , dtype=np.complex128 )
    # forst loop over offsite 
    cdef int i, j, L
    cdef double l, m, n, gsp
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
            gsp = y*pow(x,(pow(IrI,_nc)))*pow(IrI,-_n)
            # what is gsp when r=g=0?
            HGamma[4*i:4*(i+1), 4*j:4*(j+1)] = np.array([
                [ gsp*param[2],      gsp*param[3]*l,                              gsp*param[3]*m,                              gsp*param[3]*n],
                [-gsp*param[3]*l,    gsp*param[4]*l*l + gsp*param[5]*(1-l*l ), gsp*m*l*(param[4]-param[5]),                 gsp*n*l*(param[4] - param[5])],
                [-gsp*param[3]*m,    gsp*l*m*(param[4] - param[5]),               gsp*param[4]*m*m + gsp*param[5]*(1-m*m),  gsp*n*m*(param[4] - param[5])],
                [-gsp*param[3]*n,    gsp*l*n*(param[4] - param[5]),               gsp*m*n*(param[4]-param[5]),                 gsp*param[4]*n*n + gsp*param[5]*(1-n*n)]])
    if g > utils.eps:
        for i in range(L):
            l, m, n = G/g
            # what is gsp when r=g=0?
            gsp = y*pow(x,(pow(g,_nc)))*pow(g,-_n)
            HGamma[4*i:4*(i+1),4*i:4*(i+1)] = np.array([
                [ gsp*param[2],      gsp*param[3]*l,                              gsp*param[3]*m,                              gsp*param[3]*n],
                [-gsp*param[3]*l,    gsp*param[4]*l*l + gsp*param[5]*(1-l*l),  gsp*m*l*(param[4]-param[5]),                 gsp*n*l*(param[4] - param[5])],
                [-gsp*param[3]*m,    gsp*l*m*(param[4] - param[5]),               gsp*param[4]*m*m + gsp*param[5]*(1-m*m),  gsp*n*m*(param[4] - param[5])],
                [-gsp*param[3]*n,    gsp*l*n*(param[4] - param[5]),               gsp*m*n*(param[4]-param[5]),                 gsp*param[4]*n*n + gsp*param[5]*(1-n*n)]])
    else:
        for i in range(L):
            HGamma[4*i:4*(i+1),4*i:4*(i+1)] = np.diag([param[0], param[1], param[1], param[1]]) 
    return HGamma

