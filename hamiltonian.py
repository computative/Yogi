from claterkoster import HGamma
import numpy as np
from multiprocessing import Pool
import utils
import scipy as sci
import ase
from ase.lattice.cubic import Diamond

class Hamiltonian:

    def __init__(self, atoms, functions, param, gsp):

        self.N = len(atoms)*len(functions)
        a1,a2,a3 = atoms.cell[:]
        A = np.array([a1,a2,a3]).T
        
        num_neighbors = 1
        self.neighbors = []
        for i in range(-num_neighbors, num_neighbors+1):
            for j in range(-num_neighbors, num_neighbors+1):
                for k in range(-num_neighbors, num_neighbors+1):
                    self.neighbors.append(  np.dot(A, np.array([i, j, k]) ) )
        R = atoms.get_positions(wrap=True)
        
        with Pool() as pool:
            param = np.array([param[key] for key in 
                              ["Es","Ep","ss-sigma","sp-sigma","pp-sigma","pp-pi"]])
            gsp_param = np.array([gsp[key] for key in ["r0","rc","n","nc"]])
            self.HGamma = np.array(pool.starmap( HGamma, [
                (k, G,self.N, param, gsp_param, R) for k, G in enumerate(self.neighbors)] ) )

    def __call__(self, k):
        H = np.zeros( (self.N,self.N) , dtype=complex )
        for i, G in enumerate(self.neighbors):
            H += np.exp(1j*np.dot( G, k ))*self.HGamma[i]
        return H

if __name__ == "__main__":

    atoms = Diamond(
        directions=[[1,0,0], [0,1,0], [0,0,1]],
        symbol='Si',
        size=(1, 1, 1),
        latticeconstant=5.431,
        pbc=True
    )

    bowler = {  "Es": -12.2, "Ep": -5.75, 
                "ss-sigma": -1.938, "sp-sigma": 1.745,
                "pp-sigma": 3.050, "pp-pi": -1.075 }
    gsp = {"r0": 2.35, "rc": 3.8661,"n": 1.9771,"nc": 6.8702 }
    basis = ["1s","2px","2py","2pz"]
    hamiltonian = Hamiltonian(atoms, basis, bowler, gsp)

    from visualtools import VisualTools
    
    a1,a2,a3 = atoms.cell[:]
    
    sympts = VisualTools.fcc_sympts(*VisualTools.reciprocal(a1,a2,a3))
    kpath = {
        "waypts": [ sympts["K"], sympts["Gamma"],sympts["L"], sympts["K"]  ],
        "n" : [35,35,35]
    }

    eigs = []
    for i, k in enumerate(VisualTools.kpts(kpath["waypts"], kpath["n"])):
        print( "Iteration", i )
        HMatrix = hamiltonian( k )

        if np.linalg.norm(k) == 0:
            VisualTools.matrix_formatting()
            print( np.round( HMatrix,3 ))
        if not sci.linalg.ishermitian(HMatrix, atol=1e-10): 
            raise(RuntimeError("Hamiltonian nonhermitian"))
        eigs.append(np.linalg.eigvalsh( HMatrix, "L"))

    bands = np.array(eigs).T
    VisualTools.bands(bands)

