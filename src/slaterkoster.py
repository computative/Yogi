class Slaterkoster:

    def __init__(self, basis, params):
        k = len(basis.functions)
        self.prime = {}
        primes = self.prime_gen(k)
        # a map from single quantum numbers to primes
        for f_i,prime in zip(basis.functions, primes):
            self.prime[f_i] = prime

        # truncate from products of primes to a range
        self.map_prod = {}
        for i,prod in enumerate( np.unique( [i*j \
                        for i in primes for j in primes] ) ):
            self.map[prod] = i

        # truncate from single primes to a range
        self.map_single = {}
        for i, sing in enumerate( primes ) :
            self.map[sing] = i


        # sets hardcoded self._param
        from slaterkoster.param import Param

        self._param = np.zeros((k,k))

        self._param = Param.build(param)

        # sets hardcoded self._coef
        from slaterkoster.coef import Coef 
        self._coefs = Coef.build()

        # 



    # calculates
    def coef(self, d, qn1, qn2):
        dhat = d/len(d)
        l, m, n = dhat
        idx = self.map[self.prime[qn1]*self.prime[qn2]]
        return self._coefs[ idx ](l,m,n)


    # calculates
    def param(self, qn1, qn2):
        return self._params[ self.prim[qn1] , self.prim[qn2] ]


    @static
    def prime_gen(n):
        p = 3
        primes = [2]
        while len(primes) < n:
            nonprime = None
            for j in range(2,p):
                # prime if nothing divides it
                if p % j == 0:
                    nonprime = True
                    break
            if nonprime != True:
                primes.append(p)
            p += 1
        return primes

