import numpy as np
from numpy import array as ar
import utils

class Slaterkoster:
    
    def __init__(self, basis, params, neighbors):
        self.map = dict(zip( basis.functions, range(len(basis.functions)) ))
        self.L = len(self.map)
        self.neighbors = neighbors
        self.N = len(basis)
        self.K = len(self.neighbors)
        self.basis = basis
        self.param = np.array( list(params.values()) )
    #    self.F = self.build()
    #def coef(self, d, qn1, qn2):
    #    length = np.linalg.norm(d)
    #    if length == 0:
    #        if qn1 == "1s":
    #            return ar([1,0,0,0,0,0])
    #        elif qn1[:2] == "2p":
    #            return ar([0,1,0,0,0,0])
    #    l, m, n = d/length
    #    _ro = 2.35
    #    _rc = 3.8661
    #    _n = 1.9771
    #    _nc = 6.8702
    #    r = length
    #    bowler = (_ro/r)**_n*np.exp(_n* ((_ro/_rc)**_nc - (r/_rc)**_nc))
    #    if np.linalg.norm(d - [-6.78875, -9.50425, -9.50425]) < 1e-4 and qn1 == qn2 == "2px":
    #        print("test", bowler)
    #    return bowler*self.F[ self.map[qn1]][self.map[qn2] ](l,m,n)
    #def HGamma_m(self, args ):
    #    k, G, N, sk = args
    #    HGamma = np.zeros( (N,N) , dtype=complex )
    #    for i,(R_i,f_i) in enumerate(sk.basis):
    #        for j,(R_j,f_j) in enumerate(sk.basis):
    #            if (np.linalg.norm( R_j-R_i) + np.linalg.norm(G)  < utils.eps) and f_i != f_j:
    #                continue
    #            d = (R_j + G)-R_i
    #            HGamma[i,j] = np.dot(sk.coef(d, f_i,f_j ),sk.param ) 
    #    return HGamma

    def HGamma_m2(self, args ):
        k, G, N, sk, R = args
        _ro = 2.35; _rc = 3.8661; _n = 1.9771; _nc = 6.8702;
        g = np.linalg.norm(G)
        HGamma = np.zeros( (N,N) , dtype=complex )
        # forst loop over offsite 
        for i in range( len(R)):
            for j in [x for x in range(len(R)) if x != i]:
                r = R[j] + G - R[i]
                IrI = np.linalg.norm(r)
                l, m, n = r/IrI
                # what is bowler when r=g=0?
                bowler = (_ro/IrI)**_n*np.exp(_n* ((_ro/_rc)**_nc - (IrI/_rc)**_nc)) 
                HGamma[4*i:4*(i+1), 4*j:4*(j+1)] = bowler*ar([
                    [ sk.param[2],      sk.param[3]*l,                              sk.param[3]*m,                          sk.param[3]*n],
                    [-sk.param[3]*l,    sk.param[4]*l**2 + sk.param[5]*(1-l**2),    m*l*(sk.param[4]-sk.param[5]),          n*l*(sk.param[4] - sk.param[5])],
                    [-sk.param[3]*m,    l*m*(sk.param[4] - sk.param[5]),            sk.param[4]*m**2 + sk.param[5]*(1-m**2),   n*m*(sk.param[4] - sk.param[5])],
                    [-sk.param[3]*n,    l*n*(sk.param[4] - sk.param[5]),            m*n*(sk.param[4]-sk.param[5]),          sk.param[4]*n**2 + sk.param[5]*(1-n**2)]], dtype=complex)

        if g > utils.eps:
            for i in range(len(R)):
                l, m, n = G/g
                # what is bowler when r=g=0?
                bowler = (_ro/g)**_n*np.exp(_n* ((_ro/_rc)**_nc - (g/_rc)**_nc)) 
                HGamma[4*i:4*(i+1),4*i:4*(i+1)] = bowler*ar([
                    [ sk.param[2],      sk.param[3]*l,                              sk.param[3]*m,                          sk.param[3]*n],
                    [-sk.param[3]*l,    sk.param[4]*l**2 + sk.param[5]*(1-l**2),    m*l*(sk.param[4]-sk.param[5]),          n*l*(sk.param[4] - sk.param[5])],
                    [-sk.param[3]*m,    l*m*(sk.param[4] - sk.param[5]),            sk.param[4]*m**2 + sk.param[5]*(1-m**2),   n*m*(sk.param[4] - sk.param[5])],
                    [-sk.param[3]*n,    l*n*(sk.param[4] - sk.param[5]),            m*n*(sk.param[4]-sk.param[5]),          sk.param[4]*n**2 + sk.param[5]*(1-n**2)]], dtype=complex)
        else:
            for i in range(len(R)):
                HGamma[4*i:4*(i+1),4*i:4*(i+1)] = np.diag([sk.param[0], sk.param[1], sk.param[1], sk.param[1]]) 
        return HGamma




    #def build(self):
    #    global f_1s1s, f_2px2px, f_2py2py , f_2pz2pz , f_1s2px , f_2px1s, f_1s2py , f_2py1s
    #    global f_1s2pz ,f_2pz1s ,  f_2py2px , f_2px2py , f_2pz2px , f_2px2pz , f_2pz2py , f_2py2pz
    #    def f_1s1s(l,m,n):
    #        return ar([0, 0 ,  1  ,   0   ,   0   ,  0   ])
    #    def f_2px2px(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , l**2,(1-l**2) ])
    #    def f_2py2py(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , m**2,(1-m**2) ])
    #    def f_2pz2pz(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , n**2,(1-n**2) ]) 
    #    def f_1s2px(l,m,n):
    #        return ar([0, 0 ,  0  ,   l   ,   0   ,  0   ]) 
    #    def f_2px1s(l,m,n):
    #        return -ar([0, 0 ,  0  ,   l   ,   0   ,  0   ])  
    #    def f_1s2py(l,m,n):
    #        return ar([0, 0 ,  0  ,   m   ,   0   ,  0   ]) 
    #    def f_2py1s(l,m,n):
    #        return -ar([0, 0 ,  0  ,   m   ,   0   ,  0   ]) 
    #    def f_1s2pz(l,m,n):
    #        return ar([0, 0 ,  0  ,   n   ,   0   ,  0   ]) 
    #    def f_2pz1s(l,m,n):
    #        return -ar([0, 0 ,  0  ,   n   ,   0   ,  0   ])
    #    def f_2py2px(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , l*m   , -l*m ]) 
    #    def f_2px2py(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , m*l   , -m*l ]) 
    #    def f_2pz2px(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , l*n   , -l*n ])
    #    def f_2px2pz(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , n*l   , -n*l ]) 
    #    def f_2pz2py(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , m*n   , -m*n ]) 
    #    def f_2py2pz(l,m,n):
    #        return ar([0, 0 ,  0  ,   0   , n*m   , -n*m ]) 
    #    table = np.zeros((self.L,self.L)).tolist()
    #    for qn1, i in self.map.items():
    #        for qn2, j in self.map.items():
    #            table[i][j] = eval( "f_"+qn1 + qn2  )
    #    return table

