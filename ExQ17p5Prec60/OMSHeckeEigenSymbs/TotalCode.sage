## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

## load

load("Functions.sage")

## initialize

tameLevel = 68
prime = 5
G.<a,b> = DirichletGroup(68, QQ)
H = DirichletGroup(tameLevel*prime, QQ)
chi = H(a*b)
char = chi
weight = 1 - 2
eigenVal1 = 1
eigenVal2 = -1
eigenVals = [1,-1]
precPlusError = 60
a3 = 0
a7 = 0
a11 = 0
a13 = -2

##########################
## initialize for s = 1 ##
##########################

dim = 32
sign = 1

## get eigenSymbs

basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]

basis = changePrecision(basis, prec)
save(basis, 'plusBasisNewPrec')

upVecs = getUpVecs(basis)

upMat = makeUpMat(basis, upVecs)

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
save(eigenVecs, 'upPlusEigenVecs')

eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.attach(eigenSymb)
save(eigenSymbs, 'upPlusEigenSymbs')

for i in range(0, len(eigenSymbs)):
    eigenSymbs[i] = makeHeckeEigen(eigenSymbs[i], 3, a3)
save(eigenSymbs, 'heckePlusEigenSymbs')

LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = getCosetMeasures(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 5):
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_1')

LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.list()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_1')


##########################
## initialize for s = 0 ##
##########################

dim = 32
sign = -1

## get eigenSymbs

basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]
save(basis, 'minusBasisOldPrec')

basis = changePrecision(basis, prec)
save(basis, 'minusBasisNewPrec')

upVecs = getUpVecs(basis)

upMat = makeUpMat(basis, upVecs)

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
save(eigenVecs, 'upMinusEigenVecs')

eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.attach(eigenSymb)
save(eigenSymbs, 'upMinusEigenSymbs')

for i in range(0, len(eigenSymbs)):
    eigenSymbs[i] = makeHeckeEigen(phi, a3)
save(eigenSymbs, 'heckeMinusEigenSymbs')

LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = cosetIntegralsAt_s_0(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 5):
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_0')

LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.list()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_0')

