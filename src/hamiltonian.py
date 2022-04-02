from slaterkoster import sk
from numpy import np


class Hamiltonian:

    # basis  : ["1s","2s","2px","2py","2pz"]
    # params : {  "ss-sigma": 1, "sp-sigma": 1,
    #            "pp-sigma": 1, "pp-pi": 1   }

    def __init__(self, crystal, basis, params):
        self.crystal = crystal
        self.basis = basis
        self.sk = sk(basis, params)

    def __call__(self, k=0):
        n = len(self.basis)
        H = np.zeros((n,n))
        for (R_i,f_i) in self.basis:
            for (R_j,f_j) in self.basis:
                s = 0.0 + 0*1j
                for G in self.crystal.neighbors(R_j):
                    s += exp(1j*G*k)*np.dot( 
                            self.sk.coef(G-R_i, f_i,f_j ),
                            self.sk.param(f_i,f_j ) )

        return np.array(H, dtype = np.complex128)


"""

OVER HER MÅ JEG IMPLEMENTERE NERMESTE NABOAPPROKSIMASJON

"""


if __name__ == "__main__":
    from crystal import Crystal
    crl1 = Crystal(dims = (1,1,1), 
                struct = {"type":"diamond", 
                          "spp": ["Si","Si"]})

    basis = ["1s","2s","2px","2py","2pz"]
    params = {  "ss-sigma": 1, "sp-sigma": 1,
                "pp-sigma": 1, "pp-pi": 1   }

    params = {  "1s1s-sig": 1, 
         "1s2s-sig": 1, 
         "1s1p-sig": 1, 
         "1s1p-sig": 1, 
         "1s1p-sig": 1, 
         "1s1s-sig": 1, 
         "1s1s-sig": 1, 
         "1s1s-sig": 1, 
         "1s1s-sig": 1, 
       }


    H1 = Hamiltonian(crl1, basis, params)
    H1(k=0)

