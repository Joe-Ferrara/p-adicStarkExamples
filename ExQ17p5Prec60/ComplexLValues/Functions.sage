from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding

def Ideals(N, K):

## Returns a list of the ideals of norm N in a quadratic field K

    N = NN(N)
    D = K.discriminant()
    l1 = [] ## list of pairs [p dividing N, order of p dividing N]
    lprimes = N.prime_factors()
    for i in range(0, len(lprimes)):
        l1.append([lprimes[i], N.valuation(lprimes[i])])
    for i in range(0, len(lprimes)):
        if kronecker_symbol(D, l1[i][0]) == -1 and l1[i][1]%2 == 1:
            return 0
            ## if the exponent of a prime that stays prime in K is odd, then there are no primes of K with norm N
    l2 = []  ## list of ideals of K whose norm is a prime power that exactly divides n
    for i in range(0, len(l1)):
        p = l1[i][0]
        ordp = l1[i][1]
        pK = K.ideal(p)
        l3 = [] ## list of ideals above p with norm p^ordp
        if kronecker_symbol(D, p) == 0:
            scrp = pK.prime_factors()[0]
            l3.append(scrp^(ordp))
        if kronecker_symbol(D, p) == 1:
            for j in range(0, ordp + 1):
                scra = (pK.prime_factors()[0]^j)*(pK.prime_factors()[1]^(ordp - j))
                l3.append(scra)
        if kronecker_symbol(D, p) == -1:
            l3.append(pK^(ordp/2))
        l2.append(l3)
    a = l2[0]
    k = 0
    while k < len(l2) - 1:
        b = l2[k + 1]
        c = []
        for i in range(0, len(a)):
            for j in range(0, len(b)):
                c.append(a[i]*b[j])
        a = c
        k = k + 1
    return a

def Frob(K, L, P):

## K is a number field
## L is an abelian extension of K
## P is a prime of K that does not ramify
## Returns the Frobenius of L/K at P

    auts = L.automorphisms()
    PA = L.ideal(P).prime_factors()[0]
    N = P.absolute_norm()
    intbasis = L.integral_basis()
    for sig in auts:
        l1 = []
        for i in range(0, len(intbasis)):
            if (sig(intbasis[i]) - intbasis[i]^N).valuation(PA) > 0:
                l1.append(intbasis[i])
        if len(l1) == len(intbasis):
            return sig

def FindGen(K, L, b = None):

## K is a number field
## L is a cyclic extension of K
## b is a generator for the field extension L over K
## returns the first generator of Gal(L/K) that it finds

    auts = L.automorphisms()
    if b == None:
        b = L.gen()
    d = L.relative_degree()
    for sig in auts:
        l1 = []
        x = b
        for i in range(0, d-1):
            if sig(x) != b:
                l1.append(sig)
            x = sig(x)
        if len(l1) == d-1:
            return sig

def chiFrob(K, L, P, i):

## K is a field
## P is a prime of K
## L is an extension of K which is cut out by the character chi that sends a generator of Gal(L/K) found by FindGen(K, L) and sends that generator to a primitive [L:K]th root of unity.
## i is an integer, 0 leq i leq [L:K]
## Returns chi^i(Frob(P))

    d = L.relative_degree()
    F.<zeta> = CyclotomicField(d)
    D = L.relative_discriminant()
    l1 = D.prime_factors()
    for x in l1:
        if P == x:
            return 0
    gen = FindGen(K, L)
    for j in range(1, d+1):
        if Frob(K, L, P) == gen^j:
            return zeta^(i*j)

def chiArtin(K, L, A, i):

## K is a field
## A is an ideal of K
## L is an extension of K which is cut out by the character chi as in chiFrob(K, L, P, i)
## Returns chi^i(Artin(A))

    l1 = A.prime_factors()
    l2 = []
    for x in l1:
        order = A.valuation(x)
        chival = chiFrob(K, L, x, i)^(order)
        l2.append(chival)
    return prod(l2)

def coeffLfunc(n, K, L, i):

## n is a natural number
## K is a quadratic field
## L is an extension of K which is cut out by a character, say chi
## i is an integer, 0 leq i leq [L:K]
## Returns the nth coefficient of L(chi^i,s)

    if n == 1:
        return 1
    l1 = Ideals(n, K)
    if l1 == 0:
        return 0
    else:
        l2 = [chiArtin(K, L, x, i) for x in l1]
        return sum(l2)













