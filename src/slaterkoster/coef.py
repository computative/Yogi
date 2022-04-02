import numpy as np
from vector import c

class Coef:

    @static
    def build(self):

        """
        #                       1s1s-sig 1s2s-sig 1s2p-sig 2s2s-sig 2s2p-sig 2p2p-sig 2p2p-pi
        c1s1s = lambda l,m,n : c(   1   ,   0    ,   0    ,    0    ,   0   ,   0   ,   0    )
        c1s2s = lambda l,m,n : c(   0   ,   1    ,   0    ,    0    ,   0   ,   0   ,   0    )
        c2s2s = lambda l,m,n : c(   0   ,   0    ,   0    ,    1    ,   0   ,   0   ,   0    )
        c1s2x = lambda l,m,n : c(   0   ,   0    ,   l    ,    0    ,   0   ,   0   ,   0    )
        c1s2y = lambda l,m,n : c(   0   ,   0    ,   m    ,    0    ,   0   ,   0   ,   0    )
        c2s2x = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   l   ,   0   ,   0    )
        c2s2y = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   m   ,   0   ,   0    )
        c1s2z = lambda l,m,n : c(   0   ,   0    ,   n    ,    0    ,   0   ,   0   ,   0    )
        c2x2x = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   0   , l**2  ,(1-l**2))
        c2s2z = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   n   ,   0   ,   0    )
        c2x2y = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   0   , l*m   , -l*m   )
        c2y2y = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   0   , m**2  ,(1-m**2))
        c2x2z = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   0   ,  l*n  , -l*n   )
        c2y2z = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   0   ,  m*n  , - m*n  )
        c2z2z = lambda l,m,n : c(   0   ,   0    ,   0    ,    0    ,   0   ,  n**2 ,(1-n**2))

        return np.array([c1s1s, c1s2s, c2s2s, c1s2x, c1s2y, c2s2x, c2s2y, c1s2z, \
                                c2x2x, c2s2z, c2x2y, c2y2y, c2x2z, c2y2z, c2z2z])
        params = {  "1s1s-sig": 1, "1s2s-sig": 1, "1s2p-sig": 1, 
                "2s2s-sig": 1, "2s2p-sig": 1, "2p2p-sig": 1,
            "2p2p-pi": 1}        
        """

        
        """
        0 1s1s 1 1s2s 2 2s2s 3 1s2px 4 1s2py 5 2s2px 6 2s2py
        7 1s2pz 8 2px2px 9 2s2pz 10 2px2py 11 2py2py 12 2px2pz
        13 2py2pz 14 2pz2pz
        """


        #                    Es, Ep ss-sig sp-sig  pp-sig  pp-pi      
        ss = lambda l,m,n : c(0, 0 ,  1  ,   0   ,   0   ,  0   )
        sx = lambda l,m,n : c(0, 0 ,  0  ,   l   ,   0   ,  0   )
        sy = lambda l,m,n : c(0, 0 ,  0  ,   m   ,   0   ,  0   )
        sz = lambda l,m,n : c(0, 0 ,  0  ,   n   ,   0   ,  0   )
        xx = lambda l,m,n : c(0, 0 ,  0  ,   0   , l**2,(1-l**2))
        yy = lambda l,m,n : c(0, 0 ,  0  ,   0   , m**2,(1-m**2))
        zz = lambda l,m,n : c(0, 0 ,  0  ,   0   , n**2,(1-n**2))
        yx = lambda l,m,n : c(0, 0 ,  0  ,   0   , l*m   , -l*m )
        zx = lambda l,m,n : c(0, 0 ,  0  ,   0   , l*n   , -l*n )
        zy = lambda l,m,n : c(0, 0 ,  0  ,   0   , m*n   , -m*n )
        #                    s  x   y  z
        return np.array([   [ ss, sx, sy, sz],   # s
                            [-sx, xx, yx, zx],   # x
                            [-sy,-yx, yy, zy],   # y
                            [-sz,-zx,-zy, zz]]  )# z
        
        """

        1 1s
        2 2px
        3 2py
        4 2pz

        """