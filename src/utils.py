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

eps = 1e-14