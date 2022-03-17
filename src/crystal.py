import numpy as np
import matplotlib.pyplot as plt
from param import *
from linked import linked
from box import Box


T = np.array([[0.5,0.5,0.0],
           [0.0,0.5,0.5],
           [0.5,0.0,0.5]])

Tinv = np.array([[1.0,-1.0,1.0],
              [1.0,1.0,-1.0],
              [-1.0,1.0,1.0]])

sgn = lambda x: 1 if x == 0 else x/np.abs(x)
norm_T = lambda x: (4/3**.5)*np.linalg.norm(np.dot(Tinv,x))

def GSP(r,r0 = 1, rc = 2.5, n = 2, nc = 5):
    if r == 0:
        return 10^15
    return (r0/r)**2*np.exp(-n*(r/rc)**nc)*np.exp(n*(r0/rc)**nc)

# the linked datatype is the crystal datatype because 
# a crystal is set of linked attributes (coordinates)

class entries:
    def __init__(self,items):
        result = []
        for Id,links,prty in items:
            result.append({"id": Id, "links": links, "prty": prty})

# add links to nuclei            
# change coordinates of anion cation in cell
# find neighbors of the nucleus and link electrons
# make a list of all electron entries and coords


class Crystal:

    """
    Gamma pt periodic boundary condition is to set k = 0
    and solve secular eqn for 
    """

    def __init__(self,n):

        # linear transformation of a cube to
        # oblique fcc primitive cell. The
        # anion-cation bond is an eigenvector
        
        """
        Her setter vi inn endrig for å endre fasongen til 
        krystallen. Istedenfor padding, gjør sett in pbc
        """

        self.settings = {"crystal_dim": (1,1), "" }

        self.bounding_box = Box(self.settings["BoundingBox"])

        self.lattice = self._make_latice()

        add_nuclei()
        add_orbitals()

    def make_latticepts

    def add_nuclei(self):
        # list of coordinates
        attribs = []

        """
        Her må vi endre, evt fjerne molekyler 
        som er utenfor en bounding box
        """

        # insert cation and anions at each lattice vertex
        for vertex in self.lattice.reshape(self.n**3,3):
            x,y,z = vertex
            # anion coords, anion T
            attribs.append([[x,y,z],True,vertex.tolist()])
            #cation coords, anion F
            attribs.append([[x + 0.25,y + 0.25,z + 0.25],
                            False,vertex.tolist()])

        """
        Her linker vi atomer. Endre slik at atomet linker til 
        de fire nermeste atomene. Linkingen ma vaere en 
        separat metode slik at den brukes når atomer flyttes
        """

        #instanciate a linked attributes object
        self.atom = linked(attribs)
        for entry in self.atom:
        #for entry in self.atom.data:
            x, y, z = entry[2][0]
            _id = entry[0]
            if entry[2][1] == 1: # if anion
                self.atom.link_attr(_id,0,
                                    [[x+0.25,y+0.25,z+0.25],
                                     [x+0.25,y-0.25,z-0.25],
                                     [x-0.25,y-0.25,z+0.25],
                                     [x-0.25,y+0.25,z-0.25]])
            else:                # if cation
                self.atom.link_attr(_id,0,
                                    [[x-0.25,y-0.25,z-0.25],
                                     [x+0.25,y+0.25,z-0.25],
                                     [x-0.25,y+0.25,z+0.25],
                                     [x+0.25,y-0.25,z+0.25]])

    """
    Burde stå før linking. Call linkingen pa nytt etter
    at atomet er flyttet
    """

    # perturb elements in the primitive cell
    def perturb(self, cell, anion, shift):
        # get coordinate of lattice
        coord = self.lattice[cell[0]][cell[1]][cell[2]]
        #for entry in self.atom.data:
        for entry in self.atom:
            if ( all(entry[2][0] == coord) and entry[2][1] == anion):
                _id = entry[0]
        self.atom.data[_id][2][0] += np.array(shift)
        self.atom.data[_id][2][0] = self.atom.data[_id][2][0].tolist()
        
        coords = []
        links = []
        C = []
        #for entry in self.atom.data:
        for entry in self.atom:
            coord = entry[2][0]
            coords.append(coord)
            _cell = np.dot(Tinv,np.array(entry[2][2]))

            if np.all(_cell >= np.array([1,1,1])) and \
               np.all(_cell <= np.array([self.n,self.n,self.n])-2):
                C.append(True)
            else:
                C.append(False)
            _temp = []
            for _link_id in entry[1]:
                link_entry = self.atom.data[_link_id]
                link_coord = link_entry[2][0]
                xs = [coord[0], link_coord[0]]
                ys = [coord[1], link_coord[1]]
                zs = [coord[2], link_coord[2]]
                _temp.append([xs,ys,zs])
            links.append(_temp)
        coords = np.array(coords)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(coords[C,0], coords[C,1], coords[C,2], c = "b")
        for atom_links,c in zip(links,C):
            if not c:
                continue
            for link in atom_links:
                ax.plot(*link, c="b")
        plt.show()
        

    # add orbitals
    # benchmarked: - orbital count correct
    #              - linking appears ok
    #              - species appears ok
    def add_orbitals(self):
        attribs = []
        #for entry in self.atom.data:
        for entry in self.atom:
            vertex = entry[2][2]
            for state in ["1s","1px","1py","1pz","2s"]:
                attribs.append([entry[2][0], state, vertex])
        self.electron = linked(attribs)
        #for entry in self.electron.data:
        for entry in self.electron:
            # scope: find all atom links at the coordinates of entry
            atom = self.atom.get_attr(0,entry[2][0])[0]
            links = self.atom.retrieve(atom[0])
            temp = []
            for link_entry in links:
                #__r = np.array(link_entry[2][0])
                #__R = np.array(atom[2][0])
                #distance = np.linalg.norm(__r-__R)
                #if distance < 0.433012701892219 or distance > 0.43301270189222:
                #    print(distance)
                #    exit()
                temp.append(link_entry[2][0])
            # scope-end: temp holds all coords that
            #            should link to electron
            # link electron with id entry[0] to coords in temp
            self.electron.link_attr(entry[0],0,temp)
        
    # output hamiltonian
    def hamiltonian(self,wavenumber):
        orbitals = []
        coords = []
        states = []
        # loop over each axis of the crystal lattice
        # to get the state type and coords of every electron
        for i in range(1,self.n-1):
            for j in range(1,self.n-1):
                for k in range(1,self.n-1):
                    # search for all atoms associated
                    # with each primitive cells 
                    atoms = self.atom.get_attr(2,self.lattice[i,j,k])
                    # loop over the pair of atoms in the primitive cell
                    for entry in atoms:
                        # search electrons associated with each atom
                        orbital = self.electron.get_attr(0,entry[2][0])
                        states.extend(orbital)
                        for _orbital in orbital:
                            coords.append(_orbital[2][0])
                            orbitals.append(_orbital[2][1])
        self.states = np.array(orbitals)
        self.coords = np.array(coords)
        # states is a list of all electron states
        # in the interior of crystal
        n = self.n - 2
        H = []
        for _state in states:
            # set _anion to "a" if anion, else set to "c"
            _coord = _state[2][0]
            _entry = self.atom.get_attr(0,_coord)[0]
            _anion = "a" if _entry[2][1] else "c"
            # empty row
            row = []
            for state in states:
                # set anion to "a" if anion, else set to "c"
                coord = state[2][0]
                entry = self.atom.get_attr(0,coord)[0]
                anion = "a" if entry[2][1] else "c"
                string = _state[2][1] + _anion + "_" + state[2][1] + anion
                counts = [string.count("x"), \
                          string.count("y"), \
                          string.count("z")]
                trig = 0
                atom = self.atom.get_attr(0,_coord)[0]
                for __entry in self.atom.retrieve(atom[0]):
                    d = np.array(__entry[2][0]) - np.array(_state[2][0])
                    sign = 1
                    for di,xi in zip(d,counts):
                        sign *= sgn(di)**xi
                    trig += sign*np.exp(1j*np.dot(d,wavenumber))/4
                # prepare string for selecting parameter from param
                #R = norm_T(np.array(_state[2][0]) - np.array(state[2][0]))
                #print(R)
                R = np.linalg.norm(np.array(_state[2][0]) - np.array(state[2][0]))/0.4330127018922193
                #print(R)
                #print(" ")
                #if state == _state:
                if R == 0:
                    trig = 1
                    R = 1
                #print(R)
                cell = param[string]*GSP(R)*trig
                row.append(cell)
            H.append(row)
        return np.array(H, dtype = np.complex128)


