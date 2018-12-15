## this code calculates the complex LValues and checks that the ordered algebraic Stark Units are in the correct order

load("Functions.sage")

## calculating primitive L-functions

## set parameters

R.<x> = PolynomialRing(QQ)
F.<a> = NumberField(x^2 + 31)

C.<zeta3> = CyclotomicField(3)

prec = 500
CompNum = ComplexField(prec)

numberCoeffs = 10000


## make coefficients of not twisted L-functions

coeffsNotTwisted = [1]
coeffsId = [1]

## this is our choice for generator of class group
aIdeal = F.ideal(2, (1 + a)/2)

for n in range(2, numberCoeffs):
    print n
    idealsNormn = Ideals(n, F)
    if idealsNormn == 0:
        coeffsId.append(0)
        coeffsNotTwisted.append(0)
    if idealsNormn != 0:
        coef = 0
        coeffsId.append(len(idealsNormn))
        for ideal in idealsNormn:
            if ideal.is_principal() == True:
                coef = coef + 1
            if (aIdeal*ideal).is_principal() == True:
                coef = coef + zeta3^2
            if ((aIdeal^2)*ideal).is_principal() == True:
                coef = coef + zeta3
        coeffsNotTwisted.append(coef)

## make other L-function coefficients

print "making other coefficients"

coeffsPsi = []
coeffsPsi2 = []
coeffsTwistedByPsi = []
coeffsTwistedByPsi2 = []

for n in range(1, numberCoeffs):
    print n
    if n%3 == 0:
        coeffsPsi.append(0)
        coeffsPsi2.append(0)
        coeffsTwistedByPsi.append(0)
        coeffsTwistedByPsi2.append(0)
    if (n%9 == 1 or n%9 == 8):
        coeffsPsi.append(1*coeffsId[n-1])
        coeffsPsi2.append(1*coeffsId[n-1])
        coeffsTwistedByPsi.append(coeffsNotTwisted[n-1])
        coeffsTwistedByPsi2.append(coeffsNotTwisted[n-1])
    if (n%9 == 2 or n%9 == 7):
        coeffsPsi.append(zeta3*coeffsId[n-1])
        coeffsPsi2.append(zeta3^2*coeffsId[n-1])
        coeffsTwistedByPsi.append(zeta3*coeffsNotTwisted[n-1])
        coeffsTwistedByPsi2.append(zeta3^2*coeffsNotTwisted[n-1])
    if (n%9 == 4 or n%9 == 5):
        coeffsPsi.append(zeta3^2*coeffsId[n-1])
        coeffsPsi2.append(zeta3*coeffsId[n-1])
        coeffsTwistedByPsi.append(zeta3^2*coeffsNotTwisted[n-1])
        coeffsTwistedByPsi2.append(zeta3*coeffsNotTwisted[n-1])

print "ok"


## make coefficients complex numbers

allCoeffs = [coeffsId, coeffsPsi, coeffsPsi2, coeffsNotTwisted, coeffsTwistedByPsi, coeffsTwistedByPsi2]

allComplCoeffs = []
for coeffs in allCoeffs:
    complCoeffs = []
    for i in range(0, len(coeffs)):
        coeff = CompNum(coeffs[i])
        complCoeffs.append(coeff)
    allComplCoeffs.append(complCoeffs)

## make the pairs that are related by functional equation

coeffsPairs = [[allComplCoeffs[0], allComplCoeffs[0]], [allComplCoeffs[1], allComplCoeffs[2]], [allComplCoeffs[2], allComplCoeffs[1]], [allComplCoeffs[3], allComplCoeffs[3]], [allComplCoeffs[4], allComplCoeffs[5]], [allComplCoeffs[5], allComplCoeffs[4]]]

## root numbers 
rootNumbers = [1, CompNum(exp(-2*pi*I/9)), CompNum(exp(2*pi*I/9)), 1, CompNum(exp(-2*pi*I/9)), CompNum(exp(2*pi*I/9))]

conductors = [31, 31*9^2, 31*9^2, 31, 31*9^2, 31*9^2]

gammaV = [0,1]

wt = 1

poles = []

residues = 'automatic'

init = '1'

## make Lfunctions using Dokchitser

Lfunctions = []
LfuncCheckFunctEq = []
LValues_at_1 = []

## do zeta function by hand
zeta_function = F.zeta_function(prec)
Lfunctions.append(zeta_function)
LfuncCheckFunctEq.append(zeta_function.check_functional_equation())
LValues_at_1.append(CompNum(2*pi*3/(2*(31.sqrt())))) ## putting residue here


for i in range(1, 6):
    print i

    ## set parameters for Dokchitser

    cond = conductors[i]
    coeffs = coeffsPairs[i]
    root = rootNumbers[i]

    ## make the Lfunctions
    print 0
    Lfunct = Dokchitser(cond, gammaV, wt, root, poles, residues, prec, init)
    print 1
    Lfunct.init_coeffs(coeffs[0], 1, coeffs[1])
    print 2
    Lfunctions.append(Lfunct)
    print 3
    LfuncCheckFunctEq.append(Lfunct.check_functional_equation())
    print 4
    LValues_at_1.append(Lfunct(1))
    print 5

## save stuff

save(Lfunctions, "Lfunctions")
save(LfuncCheckFunctEq, "LfuncCheckFunctEq")
save(LValues_at_1, "LValues_at_1")

print "derivatives"

## make the derivatives at 0

LDer_at_0_0 = CompNum(3*((1/3).log())) ##zeta functions without Euler factor at 3

constant1 = CompNum(rootNumbers[1]*(conductors[1].sqrt())*(1/(2*pi)))
LDer_at_0_1 = constant1*LValues_at_1[2]

constant2 = CompNum(rootNumbers[2]*(conductors[2].sqrt())*(1/(2*pi)))
LDer_at_0_2 = constant2*LValues_at_1[1]

constant4 = CompNum(rootNumbers[4]*(conductors[4].sqrt())*(1/(2*pi)))
LDer_at_0_4 = constant4*LValues_at_1[5]

constant5 = CompNum(rootNumbers[5]*(conductors[5].sqrt())*(1/(2*pi)))
LDer_at_0_5 = constant5*LValues_at_1[4]

LDerivativesAt_0 = [LDer_at_0_0, LDer_at_0_1, LDer_at_0_2, 0, LDer_at_0_4, LDer_at_0_5]


## make leading terms of partial zeta values
## ordering matter

LDers0 = LDerivativesAt_0
Zeta3 = CompNum(exp(2*pi*I/3))
partialZetaVals = []
for i in range(0, 3):
    for j in range(0, 3):
        zetaVal = 0
        zetaVal = zetaVal + LDers0[0]
        zetaVal = zetaVal + 2*((Zeta3^(-i)*LDers0[1]).real_part())
        zetaVal = zetaVal + 2*(((Zeta3^(-j) + Zeta3^(-2*j))*Zeta3^(-i)*LDers0[4]).real_part())  
        zetaVal = (1/9)*zetaVal
        partialZetaVals.append(zetaVal)

orderedStarkUnitsData = load("orderedStarkUnitsData")

R.<x> = PolynomialRing(CompNum)
f = x^3 + x + 1

b = CompNum(f.roots()[0][0])
a = CompNum((-31).sqrt())
c = CompNum(exp(2*pi*I/9) + exp(-2*pi*I/9))

StUnitsNum = []
for n in range(0, 9):
    StUnitNum = CompNum(0)
    data = orderedStarkUnitsData[n]
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 2):
                StUnitNum = StUnitNum + CompNum(data[6*i + 2*j + k])*c^i*b^j*a^k
    StUnitsNum.append(StUnitNum)

StarkLogsNum = []
for i in range(0, 9):
    StarkLogsNum.append((-1)*((StUnitsNum[i].abs()).log()))

errors = []
for i in range(0, 9):
    error = (partialZetaVals[i] - StarkLogsNum[i]).abs()
    errors.append(error)





