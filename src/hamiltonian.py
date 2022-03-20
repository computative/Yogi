
class Hamiltonian:

    _instance = None

    def __new__(cls, crystal):
        if cls._instance == None:
            cls._instance = super().__new__(cls)
        cls._instance.__init__(crystal)
        return cls._instance.evaluate()

    def __init__(self, crystal):
        # initiate as normal
        print("__init__")
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
    crl1 = Crystal(dims = (1,1,1), 
                struct = {"type":"diamond", 
                          "spp": ["Si","Si"]})
    crl2 = Crystal(dims = (3,3,3), 
                struct = {"type":"diamond", 
                          "spp": ["Si","Si"]})

    H1 = Hamiltonian(crl1)
    H2 = Hamiltonian(crl2)
