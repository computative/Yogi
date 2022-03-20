import numpy as np


class Crystal:

    def __init__(self,dim,struct):

        self.settings = {"dims": dim, 
                    "struct": {"type": struct["type"], 
                               "spp": struct["spp"] },
                    "colors": {"Si": "grey"},
                    "cutoff": 1 }
        
        
        if self.settings["struct"]["type"] == "diamond":
            self.nuclei = self._diamond(self.settings["dims"], 
                                self.settings["struct"]["spp"])
        else:
            print("E: Structure not found")
    

    def _diamond(self, dims, species):

        from struc import DiamondCell as Cell
        
        # the method must contain empty dict and
        # nx,ny,nz number of blocks and start_Id

        sp1, sp2 = species
        nuclei = {}
        nuclei[sp1] = {}
        nuclei[sp2] = {}
        nx,ny,nz = dims
        start_Id = 0

        for i in range( nx ):
            for j in range( ny ):
                for k in range( nz ):
                    # end_plate sees
                    end_plate =  np.argwhere(
                        np.array([nx - i,ny - j,nz - k]) == 1 ).T[0]
                    origin = (i,j,k)
                    cell, added = Cell(origin, start_Id,
                        species, end_plate)
                    start_Id += added
                    for sp in np.unique(np.array(
                            self.settings["struct"]["spp"])):
                        nuclei[sp].update(cell[sp])
        return nuclei


    def plot(self, show_ids, show_spheres):
        # gi et bilde av krystallen med id-tekst isteden-
        # for kule. Forskjellige farger avhengig av 
        # element-type
        
        import matplotlib.pyplot as plt

        ax = plt.figure().add_subplot(projection='3d')

        for sp in np.unique(np.array(self.settings["struct"]["spp"])):
            coords = self.nuclei[sp].values()
            Ids = self.nuclei[sp].keys()
            col = self.settings["colors"][sp]
            xs = [coord[0] for coord in coords]
            ys = [coord[1] for coord in coords]
            zs = [coord[2] for coord in coords]
            for Id, x, y, z in zip(Ids, xs, ys, zs):
                label = str(Id)
                if show_ids:
                    ax.text(x, y, z, label, color = col)
                if show_spheres:
                    ax.scatter(x,y,z,color=col)

        # Tweaking display region and labels
        #ax.set_xlim(-1, 2)
        #ax.set_ylim(-1, 2)
        #ax.set_zlim(-1, 2)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        plt.show()


    def perturb(self, Id, shift):
        # Id is the position of nucleus in the list
        for sp in self.nuclei:
            for ID in self.nuclei[sp]:
                if ID == Id:
                    x = self.nuclei[sp][ID]
                    self.nuclei[sp][ID] = (
                        x[0] + shift[0],
                        x[1] + shift[1],
                        x[2] + shift[2])


if __name__ == "__main__":
    crystal = Crystal()
