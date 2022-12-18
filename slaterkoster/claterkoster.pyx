import numpy as np
from numpy import array as ar
import utils

cimport numpy as np

np.import_array()

class Slaterkoster:

    def __init__(self, basis, params, neighbors):
        self.map = dict(zip( basis.functions, range(len(basis.functions)) ))

        self.L = len(self.map)
        self.neighbors = neighbors
        self.N = len(basis)
        self.K = len(self.neighbors)
        self.basis = basis
        self.param = np.array( list(params.values()) )

        self.F = self.build()


    # calculates
    def coef(self, d, qn1, qn2):
        length = np.linalg.norm(d)
        if length == 0:
            if qn1 == "1s":
                return ar([1,0,0,0,0,0])
            elif qn1[:2] == "2p":
                return ar([0,1,0,0,0,0])
        elif length > 3.2:
            return ar([0,0,0,0,0,0])
        
        l, m, n = d/length
        
        bowler = ( 2.35/length )**1.9771
        return bowler*self.F[ self.map[qn1]][self.map[qn2] ](l,m,n)


    def HGamma(self):
        N = self.N; K = self.K
        HGamma = np.zeros( (K,N,N) , dtype=complex )
        for i,(R_i,f_i) in enumerate(self.basis):
            for j,(R_j,f_j) in enumerate(self.basis):
                for k, G in enumerate(self.neighbors):
                    if (np.linalg.norm( R_j-R_i) + np.linalg.norm(G)  < utils.eps) \
                                        and f_i != f_j:
                        continue
                    
                    HGamma[k,i,j] = np.dot(
                        self.coef(((R_j + G)-R_i), f_i,f_j ),
                        self.param ) 
                    
        return HGamma



    def build(self):

        def sgn(f):
            def minus(l,m,n):
                return -f(l,m,n)
            return minus

        #                         Es Ep ss-sig sp-sig  pp-sig  pp-pi      
        f_1s1s = lambda l,m,n : ar([0, 0 ,  1  ,   0   ,   0   ,  0   ])
        f_2px2px = lambda l,m,n : ar([0, 0 ,  0  ,   0   , l**2,(1-l**2) ])
        f_2py2py = lambda l,m,n : ar([0, 0 ,  0  ,   0   , m**2,(1-m**2) ])
        f_2pz2pz = lambda l,m,n : ar([0, 0 ,  0  ,   0   , n**2,(1-n**2) ])

        f_1s2px = lambda l,m,n : ar([0, 0 ,  0  ,   l   ,   0   ,  0   ])
        f_2px1s = sgn(f_1s2px)

        f_1s2py = lambda l,m,n : ar([0, 0 ,  0  ,   m   ,   0   ,  0   ])
        f_2py1s = sgn(f_1s2py)

        f_1s2pz = lambda l,m,n : ar([0, 0 ,  0  ,   n   ,   0   ,  0   ])
        f_2pz1s = sgn(f_1s2pz)

        f_2py2px = lambda l,m,n : ar([0, 0 ,  0  ,   0   , l*m   , -l*m ])
        f_2px2py = lambda l,m,n : ar([0, 0 ,  0  ,   0   , m*l   , -m*l ])

        f_2pz2px = lambda l,m,n : ar([0, 0 ,  0  ,   0   , l*n   , -l*n ])
        f_2px2pz = lambda l,m,n : ar([0, 0 ,  0  ,   0   , n*l   , -n*l ])

        f_2pz2py = lambda l,m,n : ar([0, 0 ,  0  ,   0   , m*n   , -m*n ])
        f_2py2pz = lambda l,m,n : ar([0, 0 ,  0  ,   0   , n*m   , -n*m ])

        table = np.zeros((self.L,self.L)).tolist()

        for qn1, i in self.map.items():
            for qn2, j in self.map.items():
                table[i][j] = eval( "f_"+qn1 + qn2  )

        #table = np.array(table)

        #table[[1,2,3],:] = table[[2,3,1],:]
        #table[:,[1,2,3]] = table[:,[2,3,1]]
        return table




class BasisIterator:
    def __init__(self, R, functions):
        cdef int n = len(R)*len(functions)
        cdef np.ndarray self.R = np.zeros( (n*n,3), dtype=np.double )
        cdef np.ndarray self.f = np.zeros( n*n, dtype=np.double )
        for i, r in enumerate(R):
            for j, f in enumerate(functions):
                self.R[i] = r 
                self.f[i] = f
        self.contents = [ (R,f)  for R in R  for f in functions ]
        cdef int self.pos = 0
        cdef int self.len = len(self.contents)

    cdef __next__(self):
        if self.pos >= self.len :
            raise StopIteration
        self.pos += 1
        return (self.R[pos], self.f[pos])


class Basis:

    def __init__(self, crystal, functions ):
        self.crystal = crystal
        self.functions = functions
        self.R = []
        for elts in crystal.nuclei.values():
            for R in elts.values():
                self.R.append(R)

    def __iter__(self):
        return BasisIterator(self.R, self.functions)

    def __len__(self):
        return len(self.crystal)*len(self.functions)