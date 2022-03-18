import numpy as np
import matplotlib.pyplot as plt
from param import *


class Crystal:

    """
    Gamma pt periodic boundary condition is to set k = 0
    and solve secular eqn for 
    """

    def __init__(self, n):

        self.settings = {"dims": (1,1,1), 
                    "struct": ["diamond","Si","Si"] }
        
        
        if self.settings["struct"][0] == "diamond":
            self.nuclei = self._diamond(dims, 
                                self.settings["struct"][1:])
        else:
            print("E: Structure not found")


    """
    Hver struktur faar en egen metode. Senere ma jeg
    forene dem i et felles rammeverk
    """

    def _diamond(self, dims, species):

        from struct import DiamondCell as Cell
        
        """the method must contain empty dict and
        nx,ny,nz number of blocks and start_Id"""

        nuclei = {}
        nx,ny,nz = dims
        start_Id = 0

        for i in range( nx ):
            for j in range( ny ):
                for k in range( nz ):
                    """ end_plate sees """
                    end_plate =  np.argwhere(
                        np.array([nx - i,ny - j,nz - k]) == 1 )
                    cell, added = Cell(start_Id,species,endplate)
                    nuclei.update(cell)
                    start_Id += added


    def plot(self):
        # gi et bilde av krystallen med id-tekst isteden-
        # for kule. Forskjellige farger avhengig av 
        # element-type
        pass
   
    # perturb elements in the primitive cell
    def perturb(self, Id, shift):
        # Id is the position of nucleus in the list
        pass

    # output hamiltonian
    def hamiltonian(self,wavenumber):
        wavenumber = 0.0
        # basis knowledge. Loop over basis elt. in 2D-loop
        # get coords of each nucleus
        for elt1 in self.basis:
            for elt2 in self.basis:
                sum = 0.0
                for G in self.nuclei:
                    sum += 
        H = []
        return np.array(H, dtype = np.complex128)


