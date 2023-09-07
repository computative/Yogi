from claterkoster import HGamma_m
import numpy as np
from numpy import array as ar
from multiprocessing import Pool
import utils

class Hamiltonian:

    def __init__(self, crystal, functions, param):
        # init

        param = np.array([param["Es"],param["Ep"],param["ss-sigma"],\
                param["sp-sigma"],param["pp-sigma"],param["pp-pi"]])
        self.crystal = crystal
        # find W by transforming cubic coords
        a1,a2,a3 = self.crystal.lattice
        A = np.array([a1,a2,a3]).T
        
        neighbors = 1
        self.neighbors=[]
        for i in range(-neighbors, neighbors+1):
            for j in range(-neighbors, neighbors+1):
                for k in range(-neighbors, neighbors+2):
                    self.neighbors.append(  np.dot(A, ar([ i, j, k]) ) )
        self.N = len(self.crystal)*len(functions)
        R = ar(list(list(self.crystal.nuclei.values())[0].values()))
        with Pool() as pool:
            self.HGamma = ar( pool.starmap( HGamma_m, [ (  k, G,self.N, param, R)  for k, G in enumerate(self.neighbors)] ) )

    def __call__(self, k):
        N = self.N
        H = np.zeros( (N,N) , dtype=complex )
        for i, G in enumerate(self.neighbors):
            H += np.exp(1j*np.dot( G, k ))*self.HGamma[i]
        return H

if __name__ == "__main__":
    from crystal import Crystal

    theta = np.pi/3
    Q = ar([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
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

    atoms = {"Si":[ np.dot(Q,ar(x)).tolist() for x in atoms["Si"] ]}
    rep = [ np.dot(Q,ar(x)).tolist() for x in rep]

    coords = {"rep" : rep, "atoms" : atoms}

    crl = Crystal(dims = (1,1,1)).from_coords(coords)

    bowler = {  "Es": -12.2, "Ep": -5.75, 
                "ss-sigma": -1.938, "sp-sigma": 1.745,
                "pp-sigma": 3.050, "pp-pi": -1.075 }


    hamiltonian = Hamiltonian(crl, ["1s","2px","2py","2pz"], bowler)

    eigs = []

    from visualtools import VisualTools
    
    a1,a2,a3 = rep[0], rep[1], rep[2]
    sympts = VisualTools.fcc_sympts(*VisualTools.reciprocal(a1,a2,a3))
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
        if not utils.isHermitian(HMatrix):
            raise(ValueError("Hamiltonian nonhermitian"))
        eigs.append(np.linalg.eigvalsh( HMatrix, "L"))

    bands = np.array(eigs).T
    VisualTools.bands(bands)

