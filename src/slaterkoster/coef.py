import numpy as np
from vector import c

class Coef:

    @staticmethod
    def build():

        def sgn(f):
            def minus(l,m,n,p):
                return -f(l,m,n,p)
            return minus

        #                         Es Ep ss-sig sp-sig  pp-sig  pp-pi      
        ss = lambda l,m,n,p : p*c(0, 0 ,  1  ,   0   ,   0   ,  0   )
        sx = lambda l,m,n,p : p*c(0, 0 ,  0  ,   l   ,   0   ,  0   )
        sy = lambda l,m,n,p : p*c(0, 0 ,  0  ,   m   ,   0   ,  0   )
        sz = lambda l,m,n,p : p*c(0, 0 ,  0  ,   n   ,   0   ,  0   )
        xx = lambda l,m,n,p : p*c(0, 0 ,  0  ,   0   , l**2,(1-l**2))
        yy = lambda l,m,n,p : p*c(0, 0 ,  0  ,   0   , m**2,(1-m**2))
        zz = lambda l,m,n,p : p*c(0, 0 ,  0  ,   0   , n**2,(1-n**2))
        yx = lambda l,m,n,p : p*c(0, 0 ,  0  ,   0   , l*m   , -l*m )
        zx = lambda l,m,n,p : p*c(0, 0 ,  0  ,   0   , l*n   , -l*n )
        zy = lambda l,m,n,p : p*c(0, 0 ,  0  ,   0   , m*n   , -m*n )
        #                  s      x   y   z
        map = np.array([[  ss   , sx, sy, sz],      # s
                        [sgn(sx), xx, yx, zx],      # x
                        [sgn(sy), yx, yy, zy],      # y
                        [sgn(sz), zx, zy, zz]]  )   # z

        map[[1,2,3],:] = map[[2,3,1],:]
        map[:,[1,2,3]] = map[:,[2,3,1]]
        return map