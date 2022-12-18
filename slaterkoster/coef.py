import numpy as np
from vector import c

class Coef:

    def build(map):

        n = len(map)

        def sgn(f):
            def minus(l,m,n):
                return -f(l,m,n)
            return minus

        #                         Es Ep ss-sig sp-sig  pp-sig  pp-pi      
        f_1s1s = lambda l,m,n : c(0, 0 ,  1  ,   0   ,   0   ,  0   )
        f_2px2px = lambda l,m,n : c(0, 0 ,  0  ,   0   , l**2,(1-l**2))
        f_2py2py = lambda l,m,n : c(0, 0 ,  0  ,   0   , m**2,(1-m**2))
        f_2pz2pz = lambda l,m,n : c(0, 0 ,  0  ,   0   , n**2,(1-n**2))

        f_1s2px = lambda l,m,n : c(0, 0 ,  0  ,   l   ,   0   ,  0   )
        f_2px1s = sgn(f_1s2px)

        f_1s2py = lambda l,m,n : c(0, 0 ,  0  ,   m   ,   0   ,  0   )
        f_2py1s = sgn(f_1s2py)

        f_1s2pz = lambda l,m,n : c(0, 0 ,  0  ,   n   ,   0   ,  0   )
        f_2pz1s = sgn(f_1s2pz)

        f_2py2px = lambda l,m,n : c(0, 0 ,  0  ,   0   , l*m   , -l*m )
        f_2px2py = lambda l,m,n : c(0, 0 ,  0  ,   0   , m*l   , -m*l )

        f_2pz2px = lambda l,m,n : c(0, 0 ,  0  ,   0   , l*n   , -l*n )
        f_2px2pz = lambda l,m,n : c(0, 0 ,  0  ,   0   , n*l   , -n*l )

        f_2pz2py = lambda l,m,n : c(0, 0 ,  0  ,   0   , m*n   , -m*n )
        f_2py2pz = lambda l,m,n : c(0, 0 ,  0  ,   0   , n*m   , -n*m )


        table = np.zeros((n,n)).tolist()

        for qn1, i in map.items():
            for qn2, j in map.items():
                table[i][j] = eval( "f_"+qn1 + qn2  )

        #table = np.array(table)

        #table[[1,2,3],:] = table[[2,3,1],:]
        #table[:,[1,2,3]] = table[:,[2,3,1]]
        return table