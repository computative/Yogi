import numpy as np
from numpy import array as ar

class Slaterkoster:

    def __init__(self, basis, params):
        k = len(basis.functions)
        self.param = np.array( list(params.values()) )
        self.map = dict(zip( basis.functions, range(len(basis.functions)) ))

        from slaterkoster.coef import Coef 
        self._coefs = Coef.build(self.map)


    # calculates
    def coef(self, d, qn1, qn2):
        length = np.linalg.norm(d)
        if length == 0:
            if qn1 == "1s":
                return ar([1,0,0,0,0,0])
            elif qn1[:2] == "2p":
                return ar([0,1,0,0,0,0])
        elif length > 3.2:
            return ar([0,0,0,0,0,0])
        
        l, m, n = d/length
        
        bowler = ( 2.35/length )**1.9771
        return bowler*self._coefs[ self.map[qn1]][self.map[qn2] ](l,m,n)





