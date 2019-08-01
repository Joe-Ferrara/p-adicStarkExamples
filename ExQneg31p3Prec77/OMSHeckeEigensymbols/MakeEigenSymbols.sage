## when running the MakeEigenSymbs.sage code, getting the hecke eigensymbol in the minus subspace with U3 eigenvalue -1 does not work.
## the problem may be that for the U3 eigensymbol we get (not an eigensymbol for the full hecke algebra) has vector of total measures with first measure 0. The goal of this code is to troubleshoot this issue

#####################################################
#####################################################
## this code produces Up eigensymbs.               ##
## it does not produce hecke eigensymbs correctly. ##
## to get hecke eigensymbs from Up eigensymbs,     ##
## you must mess around by hand to make the        ##
## Up eigensymb a hecke eigensymb                  ##
#####################################################
#####################################################


##########################################
## when running this                    ##
## after messing around for a while     ##
## got two hecke eigensymbols for       ##
## the whole hecke algebra              ##  
## when alpha = 1, the precision is 77  ##
## when alpha = -1, the precision is 75 ##
##########################################
##########################################

tameLevel = 31
prime = 3
G.<a> = DirichletGroup(31, QQ)
H = DirichletGroup(tameLevel*prime, QQ)
chi = H(a)
char = chi
weight = 1 - 2
eigenVal1 = 1
eigenVal2 = -1
eigenVals = [1,-1]
precPlusError = 80

dim = 3
sign = -1

## get eigenSymbs

print "making ordinary basis"
basisAndPrec = makeOrdBasis(tameLevel, prime, char, weight, dim, sign, precPlusError)
basis = basisAndPrec[0]
prec = basisAndPrec[1]
save(basis, 'minusBasisOldPrec')

basis = changePrecision(basis, prec)
save(basis, 'minusBasisNewPrec')

print "making upVecs"
upVecs = getUpVecs(basis)
save(upVecs, 'minsuUpVecs')

upMat = makeUpMat(basis, upVecs)
save(upMat, 'minusUpMatrix')

print "making upMinusEigenVecs"
eigenVecs = []
for i in range(0, len(eigenVals)):
    upEigenVecs = makeEigenVectors(upMat, eigenVals[i], prec, prime)
    eigenVecs.append(upEigenVecs[0])
## save(eigenVecs, 'upMinusEigenVecs')

print "making upMinusEigenSymbs"
eigenSymbs = []
for eigenVec in eigenVecs:
    eigenSymb = makeEigenSymb(eigenVec, basis)
    eigenSymbs.append(eigenSymb)
save(eigenSymbs, 'upMinusEigenSymbs')
