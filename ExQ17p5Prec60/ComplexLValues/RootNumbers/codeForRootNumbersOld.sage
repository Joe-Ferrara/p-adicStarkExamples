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
