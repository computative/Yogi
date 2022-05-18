import numpy as np
from vector import c

class Slaterkoster:

    def __init__(self, basis, params):
        k = len(basis.functions)
        self.prime = {}
        self.param = np.array( list(params.values()) )
        primes = self.prime_gen(k)
        # a map from single quantum numbers to primes
        for f_i,prime in zip(basis.functions, primes):
            self.prime[f_i] = prime

        # truncate from products of primes to a range
        self.map_prod = {}
        for i,prod in enumerate( np.unique( [i*j \
                    for i in primes for j in primes] ) ):
            self.map_prod[prod] = i

        # truncate from single primes to a range
        self.map_single = {}
        for i, sing in enumerate( primes ) :
            self.map_single[sing] = i

        # sets hardcoded self._coef
        from slaterkoster.coef import Coef 
        self._coefs = Coef.build()

    # calculates
    def coef(self, d, qn1, qn2):
        length = np.linalg.norm(d)
        if length == 0:
            if qn1 == "1s":
                return c(1,0,0,0,0,0)
            elif qn1[:2] == "2p":
                return c(0,1,0,0,0,0)
        else:
            dhat = d/length
        l, m, n = dhat
        
        #p1 = ( 2.35/length )**1.9771
        p2 = (2.35/length)**1.9771 * np.exp(1.9771* ((2.35/3.8661)**6.8702 - (length/3.8661)**6.8702))
        #print (p1, p2)
        return  self._coefs[ 
                    self.map_single[self.prime[qn1]], 
                    self.map_single[self.prime[qn2]] 
                ](l,m,n,p2)


    @staticmethod
    def prime_gen(n):
        p = 3
        primes = [2]
        while len(primes) < n:
            nonprime = None
            for j in range(2,p):
                if p % j == 0:
                    nonprime = True
                    break
            if nonprime != True:
                primes.append(p)
            p += 1
        return primes



