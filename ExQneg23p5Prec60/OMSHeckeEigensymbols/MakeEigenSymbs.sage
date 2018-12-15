## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

## I ran this and got a full eigensymbol in three of the four cases when the precision was 10

## load

## load("Functions.sage") don't need this anymore because it is loaded in master.sage

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
precPlusError = 65
a2 = -1
a3 = -1
a7 = 0
a11 = 0
a13 = -1

##########################
## initialize for s = 1 ##
##########################

dim = 5
sign = 1

## get eigenSymbs

print "start"

basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
print "madeBasis"
basis = basisAndPrec[0]
prec = basisAndPrec[1]

basis = changePrecision(basis, prec)
print "changed precision"
## save(basis, 'plusBasisNewPrec')

upVecs = getUpVecs(basis)
print "got up vecs"

upMat = makeUpMat(basis, upVecs)
print "made Up matrix"

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
## save(eigenVecs, 'upPlusEigenVecs')

eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.append(eigenSymb)
save(eigenSymbs, 'upPlusEigenSymbs')

for i in range(0, len(eigenSymbs)):
    eigenSymbs[i] = makeHeckeEigen(eigenSymbs[i], 2, a2)
save(eigenSymbs, 'heckePlusEigenSymbs')


##########################
## initialize for s = 0 ##
##########################

dim = 5
sign = -1

## get eigenSymbs

basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]
## save(basis, 'minusBasisOldPrec')

basis = changePrecision(basis, prec)
## save(basis, 'minusBasisNewPrec')

upVecs = getUpVecs(basis)

upMat = makeUpMat(basis, upVecs)

eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
## save(eigenVecs, 'upMinusEigenVecs')

eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.append(eigenSymb)
save(eigenSymbs, 'upMinusEigenSymbs')

for i in range(0, len(eigenSymbs)):
    eigenSymbs[i] = makeHeckeEigen(eigenSymbs[i], 2, a2)
save(eigenSymbs, 'heckeMinusEigenSymbs')
