###############################################
###############################################
## this code produces the OMS p-adic L-value ##
## at s = 0.                                 ##
## making the coset integrals is a very time ##
## consuming step.                           ##
###############################################
###############################################


def getExponents(p, n):

## if a is a minimal generator of (Z/p^nZ)^*, then this returns a list, where the ith element of the list is the minimal exponent e of a such that a^e is congruent to i mod p^n
## NOTE: This is a bad convention, but I can't think of a better one. If p divides i, then the ith term in the list is 0.

    gen = DirichletGroup(p^n).unit_gens()[0]
    ans = []
    for i in range(0, p^n):
        if i%p == 0:
            ans.append(0)
        else:
            for j in range(0, (p-1)*(p^(n-1))):
                if Integers(p^n)(gen^j) == Integers(p^n)(i):
                    ans.append(j)
    return ans

def getLValue(phi, cosetIntegrals, n, e):
## NOTE: This function will only work for the specific example Qsqrt17
## (Z/p^nZ)^* = Z/p-1Z x Z/p^(n-1)Z. If psi has kernel Z/p-1Z and sends the minimal generator of (Z/p^nZ)^* to the p^(n-1)st root of unity given below, then e is the exponent of psi

    p = phi.p()
    N = phi.num_moments()
    Qpadic = Qp(p, N)
    Ctilde.<zetapn> = CyclotomicField(p^n)
    g = (zetapn - 1).minpoly()
    Rp.<X> = PolynomialRing(Qpadic)
    g = Rp(g)
    Cp.<D> = Qpadic.extension(g)
    Zetapn = D + 1
    Zeta = Zetapn^p
    exponents = getExponents(p, n)
    ans = 0
    for i in range(0, p^n):
        ans = ans + (Zeta^(exponents[i]*e))*Cp(cosetIntegrals[i])
    return ans


def cosetIntegralsAt_s_0(phi, ap, n):
    ## Returns a list of the integrals \int_{i + p^nZp} 1/z dphi({0} - {infty})(z) for i in the range from 0 to p^n.
    ## The ith term in the list is 0 if p divides i.
    p = phi.p()
    M = phi.num_moments()
    ans = []
    for i in range(0, p^n):
        if i%p == 0:
            ans.append(0)
        else:
            integral = 0
            for j in range(0, M):
                integral = integral + teich(i, p, M)^(-j-1)*(-1)^j*basic_integral_new(phi, i, n, j, ap, 1)
            ans.append(integral)
    return ans

eigenVals = [1, -1]
eigenSymbs = load('minusHeckeEigenSymbs')

print "making LValues"
LValues = []
for i in range(0, len(eigenSymbs)):
    print "cosetIntegrals"
    cosetIntegrals = cosetIntegralsAt_s_0(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 5):
        print j
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_0')



LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.list()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_0')
