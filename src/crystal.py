import numpy as np
from vector import c


class Crystal:

    def __init__(self, dims):
        self.settings = {
            "spp" : [],
            "colors": {"Si": "grey"},
            "dims": dims
        }
        self.nuclei = None

    def __len__(self):
        return sum([ 1  for sp in self.settings["spp"] \
                    for coords in self.nuclei[sp] ])

    def from_struct(self, struct):

        self.settings["spp"] = \
                    list(np.unique(np.array(struct["spp"])))
        
        # setting flag that crystal was constructed from struct
        self.struct = struct

        if struct["type"] == "diamond":

            from struc.diamond import Diamond

            self.nuclei = Diamond(self.settings["dims"], 
                                struct["spp"])
        else:
            print("E: Structure not found")
        return self


    def from_coords(self, coords):
        # setting flag that crystal was constructed from struct

        self.coords = coords
        lat, atoms = coords.values()
        self.nuclei = {}
        for sp in coords["atoms"].keys():
            self.nuclei[sp] = {}
            self.settings["spp"].append(sp)
        nx,ny,nz = self.settings["dims"]

        start_Id = 0
        for i in range( nx ):
            for j in range( ny ):
                for k in range( nz ):
                    origin = i*lat[0]+j*lat[1]+k*lat[2]
                    for sp in self.settings["spp"]:
                        for coord in atoms[sp]:
                            _dict = {start_Id:coord+origin}
                            start_Id += 1
                            self.nuclei[sp].update(_dict)
        return self


    def transform(self, A):
        _coord = {}
        for sp in self.settings["spp"]:
            _coord[sp] = {}
            for Id, coords in self.nuclei[sp].items():
                _coord[sp].\
                    update({Id:np.dot(A,coords)})
        self.nuclei = _coord


    def plot(self, show_ids=False, show_spheres=True):
        # gi et bilde av krystallen med id-tekst isteden-
        # for kule. Forskjellige farger avhengig av 
        # element-type
        
        import matplotlib.pyplot as plt

        ax = plt.figure().add_subplot(projection='3d')

        for sp in np.unique(np.array(self.settings["spp"])):
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
        #ax.set_xlim(-1,2);ax.set_ylim(-1,2);ax.set_zlim(-1,2)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        plt.show()


    def perturb(self, Id, shift):
        # Id is the position of nucleus in the list
        for sp in self.nuclei:
            for ID in self.nuclei[sp]:
                if ID == Id:
                    self.nuclei[sp][ID] = \
                        self.nuclei[sp][ID] + shift



if __name__ == "__main__":

    from numpy import pi,cos,sin

    crl1 = Crystal(dims = (1,1,1)).from_struct(
                    {"type": "diamond", "spp": ["Si","Si"]})

 
    # rotation about z axis angle pi/4.

    t = pi/4

    # A in SO(3)

    A = np.array([
        [cos(t),-sin(t),0],
        [sin(t),cos(t), 0],
        [0     ,  0   , 1]
    ])

    # transform

    crl1.transform(A)
    
    crl1.plot()


    
    lat = [c(1,0,0),c(0,1,0), c(0,0,1)]
    atoms = {"Si": [c(0., 0.,0.), 1/4*c(1,1,1), 
                    c(0,0.5,0.5), c(1/4,3/4,3/4), 
                    c(0.5,0,0.5), c(3/4,1/4,3/4), 
                    c(0.5,0.5,0), c(3/4,3/4,1/4)
    ]}
    crl2 = Crystal(dims = (2,1,1)).from_coords(
                    {"lattice": lat, "atoms": atoms})

    crl2.plot()
    
