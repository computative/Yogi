
class Hamiltonian:

    _instance = None

    def __new__(cls, crystal):
        if cls._instance == None:
            # uses default new to make a demo object
            cls._instance = super().__new__(cls)
        cls._instance.__init__(crystal)
        H = cls._instance.evaluate()
        return H

    def __init__(self, crystal):
        # initiate as normal
        self.crystal = crystal

    def evaluate(self):
        for elt1 in self.basis:
            for elt2 in self.basis:
                sum = 0.0
                for G in self.nuclei:
                    if np.linalg.norm(x-y) < r:
                        sum += 1
        H = []
        return np.array(H, dtype = np.complex128)


if __name__ == "__main__":
    from crystal import Crystal
    crstl1 = Crystal(dims = (1,1,1), 
                struct = {"type":"diamond", 
                          "spp": ["Si","Si"]})
    crstl2 = Crystal(dims = (3,3,3), 
                struct = {"type":"diamond", 
                          "spp": ["Si","Si"]})

    H1 = Hamiltonian(crstl1)
    H2 = Hamiltonian(crstl2)
