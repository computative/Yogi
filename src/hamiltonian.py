from slaterkoster.slaterkoster import Slaterkoster as sk
import numpy as np
from basisIter import Basis
from vector import c
import sys
np.set_printoptions(threshold=sys.maxsize)

class Hamiltonian:

    def __init__(self, basis, params, eps = 1e-14):
        # init
        self.basis = basis
        self.crystal = self.basis.crystal
        self.sk = sk(basis, params)
        self.eps = eps
        # find W by transforming cubic coords
        a1,a2,a3 = self.crystal.lattice
        A = np.array([a1,a2,a3]).T
        
        def T(x):
            return np.dot(A,x)

        neighbors = 1
        self.neighbors=[]
        self.shift=[]

        for i in range(-neighbors, neighbors+1):
            for j in range(-neighbors, neighbors+1):
                for k in range(-neighbors, neighbors+1):
                    self.neighbors.append( T( c( i, j, k) ) )
                    self.shift.append(  c( i, j, k)  )



    def __call__(self, k=0):
        n = len(self.basis)
        H = np.zeros( (n,n), dtype=complex )
        for i,(R_i,f_i) in enumerate(self.basis):
            for j,(R_j,f_j) in enumerate(self.basis):
                s = 0 + 1j*0
                for shift, G in zip(self.shift,self.neighbors):
                    #if (i == 0 and j==8) :
                    #    print( i,j,shift,G )
                    #    print( "(R_j + G)-R_i" , (R_j + G)-R_i,"G" , G,"R_j" , R_j ,"R_i" , R_i )
                    if (np.linalg.norm( R_j-R_i) + np.linalg.norm(G)  < self.eps) \
                        and f_i != f_j:
                        d = 0.0
                    else:
                        d = np.dot( self.sk.coef(((R_j + G)-R_i), f_i,f_j ),
                                    self.sk.param  )
                    e = np.exp(1j*np.dot( G, k ))
                    #e = np.exp(2*np.pi*1j*np.dot( shift, k ))
                    s+= e*d
                    #if (i == 0 and j==8) :
                    #    print(  e*d,e, d )
                H[i,j] = s

        return np.array(H, dtype = np.complex128)

    @staticmethod
    def isHermitian(H, eps = 1e-6):
        n1 = len(H[:,0])
        n2 = len(H[0,:])
        for i in range( n1 ):
            for j in range(i+1, n2 ):
                if np.abs( H[i,j] - (H.T.conj())[i,j] ) > eps:
                    print( "Non-Hermitian at idx (i,j) = " + str( (i,j) )\
                        + "with values for H[i,j] = " + str(H[i,j]) \
                        + " and H[i,j].T.conj() = " + str((H.T.conj())[i,j]))
        return True






if __name__ == "__main__":
    from crystal import Crystal
    
    a = 5.431

    lat = [a*c(1,0,0), a*c(0,1,0), a*c(0,0,1)]

    atoms={   "Si":  [  a*c(0,0,0),
                        a*c(0.25,0.25,0.25),
                        a*c(0,0.5,0.5),
                        a*c(0.25,0.75,0.75),
                        a*c(0.5,0,0.5),
                        a*c(0.75,0.25,0.75),
                        a*c(0.5,0.5,0),
                        a*c(0.75,0.75,0.25)
            ]}


    """
    atoms = {"Si":   [a*c(0,0,0),a*c(0.5,0.5,0.5)]}
    """
    print(lat)
    print(atoms)

    crl2 = Crystal(dims = (1,1,1))
    coords = {"lat" : lat, "atoms" : atoms}
    crl2.from_coords(coords)

 

    bowler = {  "Es": -12.2, "Ep": -5.75, 
                "ss-sigma": -1.938, "sp-sigma": 1.745,
                "pp-sigma": 3.050, "pp-pi": -1.075 }


    basis = Basis(crl2, ["1s","2px","2py","2pz"])
    hamiltonian = Hamiltonian(basis, bowler_thesis)


    from plottools import PlotTools

    sympts = PlotTools.fcc_sympts(a)

    kpath = [ sympts["K"], sympts["Gamma"], sympts["L"], sympts["K"] ]
    n = [35, 35, 35]

    eigs = []
    
    testpath = 2*np.pi/a*np.array([[ 0.125, -0.375, -0.375],
                [ 0.125,  0.125,  0.125],
                [ 0.375, -0.125, -0.125],
                [ 0.375,  0.375,  0.375]])

    testpath = 2*np.pi/a*np.array([[ 0.125, -0.375, -0.375]])

    for i, k in enumerate(PlotTools.kpts(kpath, n)):
        print(np.round(np.linalg.norm(k), 2) )
        HMatrix = hamiltonian( k )
        np.set_printoptions(suppress=True)

        if np.linalg.norm(k) == 0:
            print("Hij")
            np.set_printoptions(threshold=np.inf)
            np.set_printoptions(linewidth=np.inf)
            np.set_printoptions(suppress=True)
            print( np.round( HMatrix,3 )) 
            #exit()

        if not Hamiltonian.isHermitian(HMatrix):
            raise(ValueError("Hamiltonian nonhermitian"))

        eigs.append(np.linalg.eigvalsh( HMatrix, "L"))

    bands = np.array(eigs).T
    PlotTools.bands(bands)

