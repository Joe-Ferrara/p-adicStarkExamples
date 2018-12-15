## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

## load

load("Functions.sage")

## initialize

tameLevel = 23
prime = 5
G.<a> = DirichletGroup(23, QQ)
H = DirichletGroup(tameLevel*prime, QQ)
chi = H(a)
char = chi
weight = 1 - 2
eigenVal1 = 1
eigenVal2 = -1
eigenVals = [1,-1]
precPlusError = 60
p = 5


##########################
## initialize for s = 1 ##
##########################

dim = 5
sign = 1

## get eigenSymbs

basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]
save(basis, 'plusBasisOldPrec')

basis = changePrecision(basis, prec)
save(basis, 'plusBasisNewPrec')

upVecs = getUpVecs(basis)
save(upVecs, 'plusUpVecs')

upMat = makeUpMat(basis, upVecs)
save(upMat, 'plusUpMat')

## THERE IS A PROBLEM IN THE FOLLOWING FOR EIGENVALUE -1. THE ROW REDUCING OF THE MATRIX DID NOT WORK CORRECTLY. THE CODE THOUGHT THAT THE EIGENSPACE WAS TWO DIMENSIONAL WHEN IT WAS REALLY ONLY ONE BECAUSE THE MATRIX WAS ROW REDUCED INCORRECTLY. THIS SHOULD BE FIXED. WAS ABLE TO FIND THE EIGENSYMBOL BY HAND AND SAVE IT WITH THE SAME NAME

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    temp = []
    for j in range(0, len(upEigenVecs)):
        temp.append(upEigenVecs[j])
    eigenVecs.append(temp)
save(eigenVecs, 'plusUpEigenVecs')

eigenSymbs = []
for eigenVec in eigenVecs:
    temp = []
    for coeff in eigenVec:
        eigenSymb = makeEigenSymb(coeff, basis)
        temp.append(eigenSymb)
    eigenSymbs.append(temp)
save(eigenSymbs, 'plusUpEigenSymbs')

##########################
## initialize for s = 0 ##
##########################

dim = 5
sign = -1

## get eigenSymbs

basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]
save(basis, 'minusBasisOldPrec')

basis = changePrecision(basis, prec)
save(basis, 'minusBasisNewPrec')

upVecs = getUpVecs(basis)
save(upVecs, 'minusUpVecs')

upMat = makeUpMat(basis, upVecs)
save(upMat, 'minusUpMat')

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    temp = []
    for j in range(0, len(upEigenVecs)):
        temp.append(upEigenVecs[j])
    eigenVecs.append(temp)
save(eigenVecs, 'minusUpEigenVecs')

eigenSymbs = []
for eigenVec in eigenVecs:
    temp = []
    for coeff in eigenVec:
        eigenSymb = makeEigenSymb(coeff, basis)
        temp.append(eigenSymb)
    eigenSymbs.append(temp)
save(eigenSymbs, 'minusUpEigenSymbs')