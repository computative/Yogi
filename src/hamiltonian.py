from slaterkoster.slaterkoster import Slaterkoster as sk
import numpy as np
from basisIter import Basis
from vector import c
import sys
np.set_printoptions(threshold=sys.maxsize)

class Hamiltonian:

    def __init__(self, crystal, basis, params, eps = 1e-14):
        self.crystal = crystal
        self.basis = basis
        self.sk = sk(basis, params)
        self.eps = eps

    def __call__(self, k=0):
        n = len(self.basis)
        H = np.zeros((n,n), dtype=complex )
        for i,(R_i,f_i) in enumerate(self.basis):
            for j,(R_j,f_j) in enumerate(self.basis):
                # if R_j and R_i denotes the same atom
                if np.linalg.norm(R_i - R_j) < self.eps:
                    if f_i != f_j: # set zero if excitations distinct
                        continue
                    else: # else set to a constant
                        if f_j == "1s":
                            H[i,j] = np.dot( self.sk.param,
                                        c(1,0,0,0,0,0) )
                        elif f_j[:2] == "2p":
                            H[i,j] = np.dot( self.sk.param,
                                        c(0,1,0,0,0,0) )
                        continue
    
                s = 0 + 1j*0
                neighbors = [
                    c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),
                    c(-1,0,0),c(0,-1,0),c(0,0,-1),
                    c(-1,-1,0),c(0,-1,-1),c(-1,0,-1),
                    c(1,1,0),c(0,1,1),c(1,0,1),
                    c(1,1,1),c(-1,-1,-1),
                    c(1,-1,0),c(0,1,-1),c(1,0,-1),
                    c(-1,1,0),c(0,-1,1),c(-1,0,1),
                    c(1,1,1),c(-1,-1,-1),
                    c(-1,1,1),c(1,-1,1),c(1,1,-1),
                    c(-1,-1,1),c(1,-1,-1),c(-1,1,-1)
                ]
                
                for G in neighbors:
                    e = np.exp(1j*np.dot( G, k ))
                    d = np.dot( 
                            self.sk.coef(((R_j + G)-R_i), f_i,f_j ),
                            self.sk.param 
                        )
                    s+= e*d
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




"""

OVER HER MÅ JEG IMPLEMENTERE NERMESTE NABOAPPROKSIMASJON

"""


if __name__ == "__main__":
    from crystal import Crystal
    crl1 = Crystal(dims = (1,1,1))
    crl1.from_struct(struct = {"type":"diamond", 
                          "spp": ["Si","Si"]}, noendplate=True)
 
    """
    lat = [c(0,1/2,1/2), c(1/2,0,1/2), c(1/2,1/2,0)]
    atoms = {"Si": [c(0,0,0),c(0.25,0.25,0.25)]}
    crl2 = Crystal(dims = (1,1,1))
    coords = {"lat" : lat, "atoms" : atoms}
    crl2.from_coords(coords)
    """


    #Es    Ep     hssσ  hspσ  hppσ   hppπ
    #−12.2 −5.75 −1.938 1.745 3.050 −1.075

    basis = Basis(crl1, ["1s","2px","2py","2pz"])
    # bowler Ge
    params = {  "Es": -13.88, "Ep": -6.39, 
                "ss-sigma": -1.695, "sp-sigma": 2.366,
                "pp-sigma": 2.853, "pp-pi": -0.823 }


    hamiltonian = Hamiltonian(crl1, basis, params)


    from plottools import PlotTools

    eigs = []

    # conventional unit cell induces a cubic lattice

    #    Gamma     K       L       U        X    Gamma
    kpath = [
        [0,0,0],[3*np.pi/2,3*np.pi/2,0],
        [np.pi,np.pi,np.pi],
        [np.pi/2,2*np.pi,np.pi/2],
        [0,2*np.pi,0],[0,0,0]
    ]
    n = [     48,    30,     30,      48,      48]

    for i, k in enumerate(PlotTools.kpts(kpath, n)):

        print(i)

        HMatrix = hamiltonian( k )


        if not Hamiltonian.isHermitian(HMatrix):
            raise(ValueError("Hamiltonian nonhermitian"))

        eigs.append(np.linalg.eigvalsh( HMatrix, "L"))

    bands = np.array(eigs).T
    PlotTools.bands(bands)

