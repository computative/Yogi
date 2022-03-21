
class Hamiltonian:

    #_instance = None

    #def __new__(cls, crystal, k=0):
    #    if cls._instance == None:
    #        cls._instance = super().__new__(cls)
    #    cls._instance.__init__(crystal,k)
    #    return cls._instance.evaluate()

    def __init__(self, crystal, basis):
        print("__init__")
        self.crystal = crystal
        self.basis = basis

    def __call__(self, k=0):
        n = len(self.basis)
        H = np.array(n,n)
        for (R_i,f_i) in self.basis:
            for (R_j,f_j) in self.basis:
                sum = 0.0
                for G in self.crystal.neighbors(R_i):
                    if np.linalg.norm(x-y) < r:
                        sum += 1
        return np.array(H, dtype = np.complex128)




if __name__ == "__main__":

    from crystal import Crystal
    crl1 = Crystal(dims = (1,1,1), 
                struct = {"type":"diamond", 
                          "spp": ["Si","Si"]})


    H1 = Hamiltonian(crl1, basis)
    H1(k=0)

