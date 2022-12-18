import numpy as np
from numpy import array as ar
import numpy.random as Rand



class Crystal:

    def __init__(self, dims = (1,1,1)):
        self.settings = {
            "spp" : [],
            "colors": {"Si": "grey", "Ge": "#668f8f"},
            "dims": dims
        }
        self.nuclei = None
        self.lattice = []

    def __len__(self):
        return sum([ 1  for sp in self.settings["spp"] \
                    for coords in self.nuclei[sp] ])


    def from_coords(self, coords):

        self.coords = coords
        rep, atoms = coords.values()
        rep = ar(rep)
        self.nuclei = {}
        self.lattice.extend( 
            rep* \
            np.reshape( np.repeat(np.array( 
                  self.settings["dims"] ),3) , (3,3) 
            )
        )

        for sp in coords["atoms"].keys():
            self.nuclei[sp] = {}
            self.settings["spp"].append(sp)
            atoms[sp] = ar(atoms[sp])
        nx,ny,nz = self.settings["dims"]

        start_Id = 0
        for i in range( nx ):
            for j in range( ny ):
                for k in range( nz ):
                    origin = i*rep[0]+j*rep[1]+k*rep[2]
                    for sp in self.settings["spp"]:
                        for coord in atoms[sp]:
                            _dict = {start_Id:coord+origin}
                            start_Id += 1
                            self.nuclei[sp].update(_dict)
        return self


    def from_struct(self, struct, **kwargs):

        # setting flag that crystal was constructed from struct
        self.struct = struct
        
        exec("from struc.{} import {} as StrucType".format(
            struct["type"].lower(), struct["type"].capitalize()
        ), globals())

        struc_type = StrucType(struct["spp"])
        self.from_coords({
            "rep" : struc_type.rep, 
            "atoms" : struc_type.atoms
        })
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


    # shift the elements but be sure they stay inside the bounding box

    def shift(self, std_dev):
        random = Rand.normal( 0, std_dev, size = (len(self),3) )
        for id, eps in enumerate(random):
            self.perturb(id, eps)


    def perturb(self, Id, shift):
        # Id is the position of nucleus in the list
        # m is some number that you add on to ensure that coord ends up positive
        m = 10
        for sp in self.nuclei:
            for ID in self.nuclei[sp]:
                if ID == Id:
                    coord = self.nuclei[sp][ID] + shift \
                        + m*self.lattice[0] + m*self.lattice[1] + m*self.lattice[2]
                    self.nuclei[sp][ID] = np.mod(coord,np.diag(self.lattice))
    
    def lat_const(self, axis):
        return np.linalg.norm(self.lattice[axis])



if __name__ == "__main__":


    crl1 = Crystal(dims = (1,1,1)).from_struct(
                    {"type": "diamond", "spp": ["Si","Si"]})

 
    # rotation about z axis angle pi/4.
    from numpy import pi, cos, sin

    # A in SO(3)
    t = pi/4

    A = np.array([
        [cos(t),-sin(t),0],
        [sin(t),cos(t), 0],
        [0     ,  0   , 1]
    ])

    # transform

    crl1.transform(A)
    
    crl1.plot()

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

    crl2 = Crystal(dims = (2,2,2))
    coords = {"rep" : rep, "atoms" : atoms}
    crl2.from_coords(coords)

    crl2.plot()

    print(crl1.nuclei)

