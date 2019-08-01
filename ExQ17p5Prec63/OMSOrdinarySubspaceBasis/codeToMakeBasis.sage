#######################################################
#######################################################
## this code makes a basis of the ordinary subspace. ##
## since verifying earlier calculations at s = 0,    ##
## we only make minus part of ordinary subspace.     ##
## computationally most time consuming.              ##
## when making random OMS we lose 3 digits of        ##
## of precision.                                     ##
#######################################################
#######################################################


## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

## functions

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
precPlusError = 66

## initialize for s = 0

dim = 32
sign = -1

## make the minus basis

vecsTotMeasModp = load('minusVecsTotMeasModp.sobj')
precisions = load('minusNewPrecisions.sobj')
p = prime
i = len(precisions)
while i < dim:
    print i
    pair = makeRandOrdOMS(tameLevel, prime, char, weight, sign, precPlusError)
    vec = vector(Zmod(p), pair[0].vector_of_total_measures())
    temp = copy(vecsTotMeasModp)
    temp.append(vec)
    A = Matrix(Zmod(p), vecsTotMeasModp)
    B = Matrix(Zmod(p), temp)
    if B.rank() > A.rank():
        precisions.append(pair[1])
        save(precisions, 'minusNewPrecisions')
        vecsTotMeasModp.append(vec)
        save(vecsTotMeasModp, 'minusVecsTotMeasModp')
        save(pair[0], 'phiminus' + str(i))
        i = i + 1
    del(pair)

