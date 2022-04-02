import numpy as np

"""
0 1s1s 1 1s2s 2 2s2s 3 1s2px 4 1s2py 5 2s2px 6 2s2py
7 1s2pz 8 2px2px 9 2s2pz 10 2px2py 11 2py2py 12 2px2pz
13 2py2pz 14 2pz2pz
"""    

class Coef:

    @static
    def build(self):
        #                   ss-sig  sp-sig  pp-sig  pp-pi  
        ss = lambda l,m,n : c(  1  ,   0   ,   0   ,  0   )
        sx = lambda l,m,n : c(  0  ,   l   ,   0   ,  0   )
        sy = lambda l,m,n : c(  0  ,   m   ,   0   ,  0   )
        sz = lambda l,m,n : c(  0  ,   n   ,   0   ,  0   )
        xx = lambda l,m,n : c(  0  ,   0   , l**2,(1-l**2))
        yy = lambda l,m,n : c(  0  ,   0   , m**2,(1-m**2))
        zz = lambda l,m,n : c(  0  ,   0   , n**2,(1-n**2))
        xy = lambda l,m,n : c(  0  ,   0   , l*m   , -l*m )
        xz = lambda l,m,n : c(  0  ,   0   , l*n   , -l*n )
        yz = lambda l,m,n : c(  0  ,   0   , m*n   , -m*n )

        return np.array([ss, ss, ss, sx, sy, sx, sy, sz, \
                                xx, sz, xy, yy, xz, yz, zz])

        
