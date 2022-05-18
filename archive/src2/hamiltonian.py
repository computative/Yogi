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

        self.neighbors=[c(0,0,0)]


    def __call__(self, k=0):
        n = len(self.basis)
        H = np.zeros( (n,n), dtype=complex )
        for i,(R_i,f_i) in enumerate(self.basis):
            for j,(R_j,f_j) in enumerate(self.basis):

                if np.linalg.norm(R_i - R_j) < self.eps:
                    if f_i != f_j:
                        continue
                
                s = 0 + 1j*0

                for G in self.neighbors:
                    
                    e = np.exp(1j*np.dot( G, k ))
                    
                    d = np.dot( 
                            self.sk.coef(((R_j + G)-R_i), f_i,f_j ),
                            self.sk.param
                        )
                    s+= e*complex(d)
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
    lat_dimer = [a*c(0,0.5,0.5), a*c(0.5,0,0.5), a*c(0.5,0.5,0)]
    atoms_dimer = {"Ge": [a*c(0,0,0),a*c(0.25,0.25,0.25)]}


    lat = [a*c(1,0,0), a*c(0,1,0), a*c(0,0,1)]
    atoms = {"Si": [a*c(0,0,0),a*c(0.25,0.25,0.25)]}

    print(lat)
    print(atoms)


    crl2 = Crystal(dims = (1,1,1))
    coords = {"lat" : lat, "atoms" : atoms}
    crl2.from_coords(coords)


    bowler_thesis = {   "Es": -12.2, "Ep": -5.75, 
                        "ss-sigma": -1.938, "sp-sigma": 1.745,
                        "pp-sigma": 3.050, "pp-pi": -1.075 }

    bowler = {  "Es": -13.88, "Ep": -6.39, 
                "ss-sigma": -1.695, "sp-sigma": 2.366,
                "pp-sigma": 2.853, "pp-pi": -0.823 }

    parms = {  "Es": -3.2967, "Ep": 4.6560, 
                "ss-sigma": -1.5003, "sp-sigma": 2.7986,
                "pp-sigma": 4.2541, "pp-pi": -1.6510 }

    basis = Basis(crl2, ["1s","2px","2py","2pz"])
    hamiltonian = Hamiltonian(basis, bowler_thesis)





    from plottools import PlotTools

    sympts = PlotTools.fcc_sympts(a)
    kpath = [ sympts["K"], sympts["Gamma"], sympts["L"], sympts["K"] ]
    n = [35, 35, 35]

    eigs = []
    for i, k in enumerate(PlotTools.kpts(kpath, n)):

        print(np.round(np.linalg.norm(k), 2) )

        HMatrix = hamiltonian( k )
        np.set_printoptions(suppress=True)

        if np.linalg.norm(k) == 0:
            np.set_printoptions(suppress=True)
            print( np.round( HMatrix,3 ).real) 
            #exit()

        if not Hamiltonian.isHermitian(HMatrix):
            raise(ValueError("Hamiltonian nonhermitian"))

        eigs.append(np.linalg.eigvalsh( HMatrix, "L"))

    bands = np.array(eigs).T
    PlotTools.bands(bands)

