import numpy as np
eps = 1e-14

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

def idxmax(D):
    m,n = D.shape
    i,j = (int(np.floor(np.argmax(D) / n)) , np.argmax(D) % n)
    return i,j

def approx(H, Hp, epsilon=eps):
    D = np.abs( H - Hp )
    m,n = D.shape
    if np.sum( D ) > m*n*epsilon:
        return False
    return True

def isHermitian(H, epsilon=eps):
    # is the matrix H different from its hermitian conjugate?
    if not approx( H, H.T.conj() ):
        D = np.abs( H - H.T.conj() )
        # what are the 2d argmax of the matrix D? 
        i, j = idxmax(D)
        print( "Non-Hermitian at idx (i,j) = " + str( (i,j) )\
            + " with values for H[i,j] = " + str(H[i,j]) \
            + " and H[i,j].T.conj() = " + str((H.T.conj())[i,j]))
        return False
    return True


