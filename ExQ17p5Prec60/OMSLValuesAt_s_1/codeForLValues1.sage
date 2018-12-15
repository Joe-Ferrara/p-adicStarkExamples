## note that we've implicitely chosen am embedding here for zeta25, the 25th primitive root of unity

def getCosetMeasures(phi, ap, n):
## gets the measures of the cosets of p^n in Zp_p for the measure phi(0 - infty). Returns a list, such that the ith element of the list is the measure of i + p^nZ_p. If p divides i, then the measure is taken to be zero since we restrict the measure to Z_p^*

    p = phi.p()
    ans = []
    for i in range(0, p^n):
        if i%p == 0:
            ans.append(0)
        else:
            meas = basic_integral_new(phi, i, n, 0, ap, 1)
            ans.append(meas)
    return ans

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


eigenVals = [1,-1]
eigenSymbs = load("plusHeckeEigenSymbs")

LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = getCosetMeasures(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 5):
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_1_new')

LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.list()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_1_new')