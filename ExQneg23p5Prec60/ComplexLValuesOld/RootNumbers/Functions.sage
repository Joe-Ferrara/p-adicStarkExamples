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
    gl2 = []  ## list of ideals of K whose norm is a prime power that exactly divides n
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

def Coeff(n, K, L, i):

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


def Root(K, L, e, cond, psi = None, place = None):

## K is a quadratic field
## L is an extension of K cut out by a character chi which is of mixed signature if K is real quadratic. #If d is the discriminant of K, then we make the assumption that the places of L above the embedding of K which uses the positive square root of d stay real. (I don't think that I need this.)#
## e is the exponent of chi, so 0 leq e leq [L:K]
## psi is a Dirichlet character
## cond is the conductor of the extension of K cut out by chi^e*psi
## place refers to the infinite place of K which stays real if K is real quadratic. place = 1 chooses the positive square root and place = -1 chooses the negative square root
## Returns the root number of chi^e*psi, denoted W(chi^e*psi)
## Only works when certain ideals of K are prinicpal. If v^a is a power of a prime of K that exactly divides cond, then this code works if v^a is principal. Otherwise, the code returns 0.

    if psi == None:
        G.<chi> = DirichletGroup(3, QQ)
        psi = chi^2
        H = DirichletGroup(1, QQ)
        psi = H(psi)
    f = cond
    l1 = []
    k = CyclotomicField()
    for n in range(0, len(f.prime_factors())):
        q = f.prime_factors()[n]
        a = f.valuation(q)
        l1.append(q^a)
    ##print 1
    l2 = []
    for x in l1:
        z = 1
        for y in l1:
            if y != x:
                r = y.element_1_mod(x)
                z = r*z
        l2.append(z)
    ##print 2
    Ws = []
    for n in range(0, len(l1)):
        x = l1[n]
        if x.is_principal() == False:
            return 0
        q = x.absolute_norm()
        beta = x.gens_reduced()[0]
        ##print beta
        if K.is_totally_real() == False:
            valbeta = 1
            for i in range(0, len(l1)):
                if l1[i] != x:
                    z = 0
                    for j in range(0, len(l1)):
                        if j == i:
                            z = z + beta*l2[j]
                        else:
                            z = z + l2[j]
                    valbeta = k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))*k(valbeta)
        ##print 3
        if K.is_totally_real() == True:
            d = K.discriminant()
            if Integers(2)(d) != Integers(2)(0):
                alpha1 = ((1 + d.sqrt())/2).n()
                alpha2 = ((1 - d.sqrt())/2).n()
            else:
                alpha1 = d.sqrt().n()
                alpha2 = d.sqrt().n()
            if (psi(-1) == 1 and place == 1) or (psi(-1) == -1 and place == -1):
                phi = NumberFieldEmbedding(K, RLF, alpha2)
            if (psi(-1) == -1 and place == 1) or (psi(-1) == 1 and place == -1):
                phi = NumberFieldEmbedding(K, RLF, alpha1)
            if phi(beta) > 0:
                valbeta = 1
            if phi(beta) < 0:
                valbeta = -1
            for i in range(0, len(l1)):
                if l1[i] != x:
                    z = 0
                    for j in range(0, len(l1)):
                        if j == i:
                            z = z + beta*l2[j]
                        else:
                            z = z + l2[j]
                    if phi(z) > 0:
                        valbeta = k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))*k(valbeta)
                    if phi(z) < 0:
                        valbeta = -1*k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))*k(valbeta)
        ##print valbeta
        l3 = []
        l4 = []
        v = f.prime_factors()[n]
        if v.absolute_norm().is_prime() == False:
            p2 = v.absolute_norm()
            a = f.valuation(v)
            invres = list(x.invertible_residues())
            for y in invres:
                z = 0
                for i in range(0, len(l1)):
                    if l1[i] == x:
                        z = z + y*l2[i]
                    else:
                        z = z + l2[i]
                if K.is_totally_real() == False:
                    val = k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))
                if K.is_totally_real() == True:
                    if phi(z) > 0:
                        val = k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))
                    if phi(z) < 0:
                        val = -1*k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))
                l3.append(val)
                w = (((beta.norm())/(beta.norm().abs()))*y*(beta.galois_conjugate())).trace()
                l4.append((k.gen(p2^a))^w)
        ##print 5
        if v.absolute_norm().is_prime() == True:
            p = v.absolute_norm()
            a = f.valuation(v)
            invres = []
            for i in range(1, NN(x.absolute_norm())):
                if Integers(p)(i) != Integers(p)(0):
                    invres.append(K(i))
            for y in invres:
                z = 0
                for i in range(0, len(l1)):
                    if l1[i] == x:
                        z = z + y*l2[i]
                    else:
                        z = z + l2[i]
                if K.is_totally_real() == False:
                    val = k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))
                if K.is_totally_real() == True:
                    if phi(z) > 0:
                        val = k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))
                    if phi(z) < 0:
                        val = -1*k(psi(K.ideal(z).absolute_norm()))*k(chiArtin(K, L, K.ideal(z), e))
                l3.append(val)
                w = y
                l4.append((k.gen(p^a))^w)
        ##print 6
        l5 = []
        for i in range(0, len(l3)):
            l5.append(l3[i]*l4[i])
        ans = (valbeta/q.sqrt())*sum(l5)
        Ws.append(ans)
        ##print CC(ans)
    if K.is_totally_real() == True:
        Ws.append(k.gen(4)^3)
    D = K.absolute_different()
    WD = k(psi(D.absolute_norm()))*k(chiArtin(K, L, D, e))
    Ws.append(WD)
    return prod(Ws)













