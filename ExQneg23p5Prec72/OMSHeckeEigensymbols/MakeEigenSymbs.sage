## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

## I ran this and got a full eigensymbol in three of the four cases when the precision was 10

## load

## load("Functions.sage") don't need this anymore because it is loaded in master.sage

## initialize

#####################################################
#####################################################
## this code produces Up eigensymbs.               ##
## it does not produce hecke eigensymbs correctly. ##
## to get hecke eigensymbs from Up eigensymbs,     ##
## you must mess around by hand to make the        ##
## Up eigensymb a hecke eigensymb                  ##
#####################################################
#####################################################


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
precPlusError = 75
a2 = -1
a3 = -1
a7 = 0
a11 = 0
a13 = -1

##########################
## initialize for s = 1 ##
###############################
## for paper only need s = 0 ##
##   deleted s = 1 code      ##
##   still in lower precision##
##   folder                  ##
###############################


##########################
## initialize for s = 0 ##
##########################

dim = 5
sign = -1

## get eigenSymbs

print "making basis"
basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]
## save(basis, 'minusBasisOldPrec')

print "changing precision"
basis = changePrecision(basis, prec)
## save(basis, 'minusBasisNewPrec')

print "getting upVecs"
upVecs = getUpVecs(basis)

print "making upMatrix"
upMat = makeUpMat(basis, upVecs)

print "getting up eigenvecs"
eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
## save(eigenVecs, 'upMinusEigenVecs')

print "making up eigenvecs into eigensymbols"
eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.append(eigenSymb)
save(eigenSymbs, 'upMinusEigenSymbs')

for i in range(0, len(eigenSymbs)):
    eigenSymbs[i] = makeHeckeEigen(eigenSymbs[i], 2, a2)
save(eigenSymbs, 'heckeMinusEigenSymbs')
