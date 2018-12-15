## BEFORE RUNNING MUST LOAD master.sage FROM DIFFERENT DIRECTORY

## functions

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

## code

basis = load('plusBasisNewPrec')

upMat = load('plusUpMat')

eigenVals = [1,-1]
prec = 57
prime = 5

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs)
save(eigenVecs, 'upPlusEigenVecsNew')

eigenSymbs = []
for eigenVec in eigenVecs:
    for coeffs in eigenVec:
        eigenSymb = makeEigenSymb(coeffs, basis)
        eigenSymbs.append(eigenSymb)
save(eigenSymbs, 'upPlusEigenSymbsNew')