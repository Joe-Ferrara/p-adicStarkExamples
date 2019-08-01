#####################################################
#####################################################
## this code produces Up eigensymbs.               ##
## it does not produce hecke eigensymbs correctly. ##
## to get hecke eigensymbs from Up eigensymbs,     ##
## you must mess around by hand to make the        ##
## Up eigensymb a hecke eigensymb                  ##
#####################################################
#####################################################


def getUpVecs(basis):
    upVecs = []
    p = basis[0].p()
    for i in range(0, len(basis)):
        phi = basis[i].hecke(p)
        upVecs.append(phi)
    return upVecs


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

def columnReduce(A, p):
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

precisions = load('minusNewPrecisions')
prec = min(precisions)
print "the precision is"
print prec

basis = []
for i in range(0, 32):
    phi = load('phiminus' + str(i))
    phi = phi.change_precision(prec)
    basis.append(phi)

A = Matrix(Zmod(5), [phi.vector_of_total_measures() for phi in basis])
print "the following should be 32"
print A.rank()

save(basis, 'minusBasisNewPrec')

p = 5
prime = 5
eigenVals = [1, -1]

print "getting upVecs"
upVecs = getUpVecs(basis)
save(upVecs, 'minusUpVecs')

print "making upMat"
upMat = makeUpMat(basis, upVecs)
save(upMat, 'minusUpMat')

print "making eigenvectors"
eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
    eigenVecs.append(upEigenVecs[1])
save(eigenVecs, 'upMinusEigenVecs')

print "making eigenSymbs"
eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.append(eigenSymb)
save(eigenSymbs, 'upMinusEigenSymbs')

print "checking eigenSymbs"
check = []
for phi in eigenSymbs:
    check.append(phi.is_Tq_eigen(p))
save(check, 'check')
