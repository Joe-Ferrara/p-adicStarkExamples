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
precPlusError = 60

## initialize for s = 0

dim = 32
sign = -1

## get eigenSymbs

i = 0
vecsTotMeasModp = []
precisions = []
p = prime
while i < dim:
    pair = makeRandOrdOMS(tameLevel, prime, char, weight, sign, precPlusError)
    phi = pair[0]
    vec = vector(Zmod(p), phi.vector_of_total_measures())
    prec = pair[1]
    temp = copy(vecsTotMeasModp)
    temp.append(vec)
    A = Matrix(Zmod(p), vecsTotMeasModp)
    B = Matrix(Zmod(p), temp)
    if B.rank() > A.rank():
        precisions.append(prec)
        save(precisions, 'newPrecisions')
        vecsTotMeasModp.append(vec)
        save(phi, 'phi' + str(i))
        i = i + 1