from slaterkoster.slaterkoster import Slaterkoster as sk
import numpy as np
from slaterkoster.basisIter import Basis
from numpy import array as ar
from multiprocessing import Pool
import utils



class Hamiltonian:

    def __init__(self, crystal, functions, params):
        # init
        self.crystal = crystal

        # find W by transforming cubic coords
        a1,a2,a3 = self.crystal.lattice
        A = np.array([a1,a2,a3]).T

        neighbors = 1
        
        self.neighbors=[]

        for i in range(-neighbors, neighbors+1):
            for j in range(-neighbors, neighbors+1):
                for k in range(-neighbors, neighbors+1):
                    self.neighbors.append(  np.dot(A, ar([ i, j, k]) ) )

        self.basis = Basis(crystal, functions)
        self.sk = sk(self.basis, params, self.neighbors)
        self.N = len(self.basis)
        with Pool() as pool:
            self.HGamma = ar( pool.map( self.sk.HGamma_m, [ (  k, G,self.N, self.sk)  for k, G in enumerate(self.sk.neighbors)] ) )


    def __call__(self, k):
        N = self.N
        H = np.zeros( (N,N) , dtype=complex )
        for i, G in enumerate(self.neighbors):
            H += np.exp(1j*np.dot( G, k ))*self.HGamma[i]
        return H

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
    rep = [[a,0,0], [0,a,0], [0,0,a]]
    atoms = {"Si": [ 
        [0,0,0],
        [a/4,a/4,a/4],
        [0,a/2,a/2],
        [a/4,3*a/4,3*a/4],
        [a/2,0,a/2],
        [3*a/4,a/4,3*a/4],
        [a/2,a/2,0],
        [3*a/4,3*a/4,a/4]
    ]}
    coords = {"rep" : rep, "atoms" : atoms}

    crl = Crystal(dims = (1,1,1)).from_coords(coords)

    bowler = {  "Es": -12.2, "Ep": -5.75, 
                "ss-sigma": -1.938, "sp-sigma": 1.745,
                "pp-sigma": 3.050, "pp-pi": -1.075 }


    hamiltonian = Hamiltonian(crl, ["1s","2px","2py","2pz"], bowler)

    eigs = []



    from visualtools import VisualTools
    sympts = VisualTools.fcc_sympts(a)
    kpath = {
        "waypts": [ sympts["K"], sympts["Gamma"],sympts["L"], sympts["K"]  ],
        "n" : [35,35,35]
    }


    for i, k in enumerate(VisualTools.kpts(kpath["waypts"], kpath["n"])):
        print( "Iteration", i )
        HMatrix = hamiltonian( k )

        if np.linalg.norm(k) == 0:
            VisualTools.matrix_formatting()
            print( np.round( HMatrix,3 ))
        if not Hamiltonian.isHermitian(HMatrix):
            raise(ValueError("Hamiltonian nonhermitian"))
        eigs.append(np.linalg.eigvalsh( HMatrix, "L"))

    bands = np.array(eigs).T
    VisualTools.bands(bands)

