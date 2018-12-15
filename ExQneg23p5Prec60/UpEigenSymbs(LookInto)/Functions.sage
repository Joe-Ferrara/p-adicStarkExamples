def plus_minus_projection(phi, sign):
    if (sign == 1):
        return phi.plus_part()
    if (sign == -1):
        return phi.minus_part()

def makeRandOrdOMS(tameLevel, prime, char, weight, sign, precPlusError):
    phi = random_OMS(tameLevel*prime, prime, weight, precPlusError, char)
    phi = phi.project_to_ordinary_subspace()
    phi = plus_minus_projection(phi, sign)
    val = phi.valuation()
    phi = phi.scale(prime**(-phi.valuation()))
    newPrec = precPlusError - val
    return [phi, newPrec]

def addPhi(basis, phi):
    test = [psi for psi in basis]
    test.append(phi)
    p = phi.p()
    basisMat = Matrix(Zmod(p), [psi.vector_of_total_measures() for psi in basis])
    testMat = Matrix(Zmod(p), [psi.vector_of_total_measures() for psi in test])
    if (testMat.rank() > basisMat.rank()):
        return True
    else:
        return False

def maxList(listNums):
    n = len(listNums)
    ans = listNums[0]
    for i in range(1, n):
        if listNums[i] > ans:
            ans - listNums[i]
    return ans

def makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError):
    basis = []
    precisions = []
    while(len(basis)) < dim:
        phiAndPrec = makeRandOrdOMS(tameLevel, prime, char, weight, sign, precPlusError)
        phi = phiAndPrec[0]
        prec = phiAndPrec[1]
        if (addPhi(basis, phi) == True):
            basis.append(phi)
            precisions.append(prec)
    prec = maxList(precisions)
    return [basis, prec]

def getUpVecs(basis):
    upVecs = []
    p = basis[0].p()
    for i in range(0, len(basis)):
        phi = basis[i].hecke(p)
        upVecs.append(phi)
    return upVecs

def changePrecision(listOMS, prec):
    ans = [phi.change_precision(prec) for phi in listOMS]
    return ans

def makeUpMat(basis, upVecs):
    p = basis[0].p()
    prec = basis[0].num_moments()
    basisMat = Matrix(Zmod(p**prec), [phi.vector_of_total_measures() for phi in basis])
    matRows = []
    for n in range(0, len(basis)):
        b = vector(Zmod(p**(prec)), upVecs[n].vector_of_total_measures())
        x = basisMat.solve_left(b)
        matRows.append(x)
    ans = Matrix(Zmod(p**prec), [v for v in matRows])
    return ans

def columnReduce(A, p): ## THERE IS A PROBLEM WITH THISe
    n = A.dimensions()[0]
    nonIdRows = []
    for i in range(0, n):
        j = i
        while j < (n - 1) and Zmod(p)(A.row(i)[j]) == 0:
            j = j + 1
        if j == (n - 1) and Zmod(p)(A.row(i)[j]) == 0:
            nonIdRows.append(i)
        else:
            A.swap_columns(i, j)
            A.set_col_to_multiple_of_col(i, i, A.row(i)[i]**(-1))
            for k in range(0, i):
                A.add_multiple_of_column(k, i, -A.row(i)[k])
            for k in range(i + 1, n):
                A.add_multiple_of_column(k, i, -A.row(i)[k])
    return [A,nonIdRows]

def makeEigenVectors(mat, eigenVal, N, p): ## given a square matrix, Mat in Zmod(p^N) and an eigenvalue eigenVal in ZZ, find a basis of eigenvectors for the eigenspace with eigenvalue a. The eigenvectors are returned as lists, not vectors
    n = mat.dimensions()[0]
    identity = Matrix(Zmod(p**N), matrix.identity(n))
    A = mat - eigenVal*identity
    A = columnReduce(A, p)[0]
    nonIdRows = columnReduce(A, p)[1]
    eigenVectors = []
    for i in range(0, len(nonIdRows)):
        eigenVector = []
        for k in range(0, nonIdRows[0]):
            eigenVector.append(-A.row(nonIdRows[i])[k])
        if nonIdRows[i] == nonIdRows[0]:
            eigenVector.append(1)
        if nonIdRows[i] != nonIdRows[0]:
            eigenVector.append(0)
        for j in range(1, len(nonIdRows)):
            for k in range(nonIdRows[j-1] + 1, nonIdRows[j]):
                eigenVector.append(-A.row(nonIdRows[i])[k])
            if nonIdRows[i] == nonIdRows[j]:
                eigenVector.append(1)
            if nonIdRows[i] != nonIdRows[j]:
                eigenVector.append(0)
        if nonIdRows[-1] < (n - 1):
            for k in range(nonIdRows[-1] + 1, n):
                eigenVector.append(-A.row(nonIdRows[i])[k])
        eigenVectors.append(eigenVector)
    return eigenVectors

def makeEigenSymb(coeffs, basis):
    coeffs = [ZZ(x) for x in coeffs]
    eta = basis[0].scale(coeffs[0])
    for i in range(1, len(coeffs)):
        eta = eta + basis[i].scale(coeffs[i])
    return eta

def makeHeckeEigen(phi, q, aq):
## phi is a Up eigensymbol. This returns a Tq eigensymbol with eigenvalue aq
## NOTE: the Up eigenspace is assumed to be 2-dimensional here
## NOTE: THIS ONLY WORKS FOR THIS SPECIFIC EXAMPLES
    if phi.is_Tq_eigen(q)[0] == True:
        return phi
    else:
        phi1 = phi.hecke(q)
        phi2 = phi1.hecke(q)
        phi3 = phi2.hecke(q)
        psi = phi.scale(aq) + phi1.scale(-1)
        psi1 = phi1.scale(aq) + phi2.scale(-1)
        psi2 = phi2.scale(aq) + phi3.scale(-1)
        N = phi.num_moments()
        p = phi.p()
        b = vector(Zmod(p^N), psi2.vector_of_total_measures())
        A = Matrix(Zmod(p^N), [psi.vector_of_total_measures(), psi1.vector_of_total_measures()])
        x = A.solve_left(b)
        c0 = ZZ(x[0])
        c1 = ZZ(x[1])
        ans = phi2 + phi1.scale(-c1) + phi.scale(-c0)
        return ans


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


def cosetIntegralsAt_s_0(phi, ap, n):
    p = phi.p()
    M = phi.num_moments()
    ans = []
    for i in range(0, p^n):
        if i%p == 0:
            ans.append(0)
        else:
            integral = 0
            for j in range(0, M):
                integral = integral + teich(i, p, M)^(-j-1)*(-1)^j*basic_integral_new(phi, i, n, j, 1, 1)
            ans.append(integral)
    return ans


