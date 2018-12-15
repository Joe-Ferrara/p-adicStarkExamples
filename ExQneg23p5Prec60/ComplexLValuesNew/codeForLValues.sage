load("Functions.sage")

## calculating primitive L-functions

## set parameters

R.<x> = PolynomialRing(QQ)
F.<a> = NumberField(x^2 + 23)

C.<zeta15> = CyclotomicField(15)
zeta3 = zeta15^5
zeta5 = zeta15^3

prec = 500
CompNum = ComplexField(prec)

numberCoeffs = 10000


## make coefficients of not twisted L-functions

coeffsNotTwisted = [1]
coeffsId = [1]

## this is our choice for generator of class group
aIdeal = F.ideal(2, (a-1)/2)

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

coeffsPsis = []
coeffsTwistPsis = []

exponents = []
for n in range(0, 25):
    if n%5 == 0:
        exponents.append(0)
    else:
        for i in range(0, 20):
            if (2^i)%25 == n:
                exponents.append(i)

for j in range(1, 5):
    coeffsPsij = []
    coeffsTwistPsij = []
    for n in range(1, numberCoeffs):
        if n%5 == 0:
            coeffsPsij.append(0)
            coeffsTwistPsij.append(0)
        else:
            e = exponents[n%25]
            coeffsPsij.append(zeta5^(j*e)*coeffsId[n-1])
            coeffsTwistPsij.append(zeta5^(j*e)*coeffsNotTwisted[n-1])
    coeffsPsis.append(coeffsPsij)
    coeffsTwistPsis.append(coeffsTwistPsij)


print "ok"


## make coefficients complex numbers

allCoeffs = [coeffsId]
for i in range(0, 4):
    allCoeffs.append(coeffsPsis[i])

allCoeffs.append(coeffsNotTwisted)

for i in range(0, 4):
    allCoeffs.append(coeffsTwistPsis[i])


allComplCoeffs = []
for coeffs in allCoeffs:
    complCoeffs = []
    for i in range(0, len(coeffs)):
        coeff = CompNum(coeffs[i])
        complCoeffs.append(coeff)
    allComplCoeffs.append(complCoeffs)

## make the pairs that are related by functional equation

coeffsPairs = [[allComplCoeffs[0], allComplCoeffs[0]], [allComplCoeffs[1], allComplCoeffs[4]], [allComplCoeffs[4], allComplCoeffs[1]], [allComplCoeffs[2], allComplCoeffs[3]], [allComplCoeffs[3], allComplCoeffs[2]], [allComplCoeffs[5], allComplCoeffs[5]], [allComplCoeffs[6], allComplCoeffs[9]], [allComplCoeffs[9], allComplCoeffs[6]], [allComplCoeffs[7], allComplCoeffs[8]], [allComplCoeffs[8], allComplCoeffs[7]]]

## root numbers 

w1 = CompNum(-0.42577929156507264886250244574425170397997304183255838987704887396972426008015608391913146021583205245804556526552595874015009085520471980717758846330765573885456029513438611331236805611320626806933171 + 0.90482705246601952771366864793269759397041391104248708676737128360244744091248234119457139508611320860990466464374264360344986263141635507701405858291955566371062246154893175269173128469294675447350562*I)

w2 = CompNum(0.062790519529313376076178224565631133122484831906666272763040759631813526760782956677445714250218053998233833042175559896896940845347908639905775803333300837498542966063602984799187183853447143567481467 + 0.99802672842827156195233680686345055333690521970154311953512440681488813141422662367509455701853766530310506124950373678404481687460493653687984703827696053399723901532994570365913993041744938852675576*I)

w3 = CompNum(0.062790519529313376076178224565631133122484831906666272763040759631813526760782956677445714250218053998233833042175559896896940845347908639905775803333300837498542966063602984799187183853447143567481467 - 0.99802672842827156195233680686345055333690521970154311953512440681488813141422662367509455701853766530310506124950373678404481687460493653687984703827696053399723901532994570365913993041744938852675576*I)

w4 = CompNum(-0.42577929156507264886250244574425170397997304183255838987704887396972426008015608391913146021583205245804556526552595874015009085520471980717758846330765573885456029513438611331236805611320626806933171 - 0.90482705246601952771366864793269759397041391104248708676737128360244744091248234119457139508611320860990466464374264360344986263141635507701405858291955566371062246154893175269173128469294675447350562*I)

w5 = CompNum(1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000)

w6 = CompNum(-0.42577929156507264886250244574425170397997304183255838987704887396972426008015608391913146021583205245804556526552595874015009085520471980717758846330765573885456029513438611331236805611320626806933171 + 0.90482705246601952771366864793269759397041391104248708676737128360244744091248234119457139508611320860990466464374264360344986263141635507701405858291955566371062246154893175269173128469294675447350562*I)

w7 = CompNum(0.062790519529313376076178224565631133122484831906666272763040759631813526760782956677445714250218053998233833042175559896896940845347908639905775803333300837498542966063602984799187183853447143567481467 + 0.99802672842827156195233680686345055333690521970154311953512440681488813141422662367509455701853766530310506124950373678404481687460493653687984703827696053399723901532994570365913993041744938852675576*I)

w8 = CompNum(0.062790519529313376076178224565631133122484831906666272763040759631813526760782956677445714250218053998233833042175559896896940845347908639905775803333300837498542966063602984799187183853447143567481467 - 0.99802672842827156195233680686345055333690521970154311953512440681488813141422662367509455701853766530310506124950373678404481687460493653687984703827696053399723901532994570365913993041744938852675576*I)

w9 = CompNum(-0.42577929156507264886250244574425170397997304183255838987704887396972426008015608391913146021583205245804556526552595874015009085520471980717758846330765573885456029513438611331236805611320626806933171 - 0.90482705246601952771366864793269759397041391104248708676737128360244744091248234119457139508611320860990466464374264360344986263141635507701405858291955566371062246154893175269173128469294675447350562*I)

w10 = CompNum(1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000)

w11 = CompNum(-0.42577929156507264886250244574425170397997304183255838987704887396972426008015608391913146021583205245804556526552595874015009085520471980717758846330765573885456029513438611331236805611320626806933171 + 0.90482705246601952771366864793269759397041391104248708676737128360244744091248234119457139508611320860990466464374264360344986263141635507701405858291955566371062246154893175269173128469294675447350562*I)

w12 = CompNum(0.062790519529313376076178224565631133122484831906666272763040759631813526760782956677445714250218053998233833042175559896896940845347908639905775803333300837498542966063602984799187183853447143567481467 + 0.99802672842827156195233680686345055333690521970154311953512440681488813141422662367509455701853766530310506124950373678404481687460493653687984703827696053399723901532994570365913993041744938852675576*I)

w13 = CompNum(0.062790519529313376076178224565631133122484831906666272763040759631813526760782956677445714250218053998233833042175559896896940845347908639905775803333300837498542966063602984799187183853447143567481467 - 0.99802672842827156195233680686345055333690521970154311953512440681488813141422662367509455701853766530310506124950373678404481687460493653687984703827696053399723901532994570365913993041744938852675576*I)

w14 = CompNum(-0.42577929156507264886250244574425170397997304183255838987704887396972426008015608391913146021583205245804556526552595874015009085520471980717758846330765573885456029513438611331236805611320626806933171 - 0.90482705246601952771366864793269759397041391104248708676737128360244744091248234119457139508611320860990466464374264360344986263141635507701405858291955566371062246154893175269173128469294675447350562*I)

w15 = CompNum(1)

rootNumbers = [w15, w3, w12, w6, w9, w5, w8, w2, w11, w14]

conductors = [23, 23*25^2, 23*25^2, 23*25^2, 23*25^2, 23, 23*25^2, 23*25^2, 23*25^2, 23*25^2]

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
LValues_at_1.append(CompNum(2*pi*3/(2*(23.sqrt())))) ## putting residue here


for i in range(1, 10):
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

LDer_at_0_0 = CompNum(3*((1/5).log())) ##zeta functions without Euler factor at 3

constant1 = CompNum(rootNumbers[1]*(conductors[1].sqrt())*(1/(2*pi)))
LDer_at_0_1 = constant1*LValues_at_1[2]

constant2 = CompNum(rootNumbers[2]*(conductors[2].sqrt())*(1/(2*pi)))
LDer_at_0_2 = constant2*LValues_at_1[1]

constant3 = CompNum(rootNumbers[3]*(conductors[3].sqrt())*(1/(2*pi)))
LDer_at_0_3 = constant3*LValues_at_1[4]

constant4 = CompNum(rootNumbers[4]*(conductors[4].sqrt())*(1/(2*pi)))
LDer_at_0_4 = constant4*LValues_at_1[3]

constant5 = CompNum(rootNumbers[5]*(conductors[5].sqrt())*(1/(2*pi)))
LDer_at_0_5 = constant5*LValues_at_1[5]

constant6 = CompNum(rootNumbers[6]*(conductors[6].sqrt())*(1/(2*pi)))
LDer_at_0_6 = constant6*LValues_at_1[7]

constant7 = CompNum(rootNumbers[7]*(conductors[7].sqrt())*(1/(2*pi)))
LDer_at_0_7 = constant7*LValues_at_1[6]

constant8 = CompNum(rootNumbers[8]*(conductors[8].sqrt())*(1/(2*pi)))
LDer_at_0_8 = constant8*LValues_at_1[9]

constant9 = CompNum(rootNumbers[9]*(conductors[9].sqrt())*(1/(2*pi)))
LDer_at_0_9 = constant9*LValues_at_1[8]





LDerivativesAt_0 = [LDer_at_0_0, LDer_at_0_1, LDer_at_0_2, LDer_at_0_3, LDer_at_0_4, 0, LDer_at_0_6, LDer_at_0_7, LDer_at_0_8, LDer_at_0_9]

## reorder so in the order of powers of chipsi

OrderedLDersAt_0 = [LDer_at_0_0, LDer_at_0_6, LDer_at_0_8, LDer_at_0_4, LDer_at_0_7, 0, LDer_at_0_1, LDer_at_0_8, LDer_at_0_9, LDer_at_0_2, 0, LDer_at_0_6, LDer_at_0_3, LDer_at_0_9, LDer_at_0_7]


## make leading terms of partial zeta values
## ordering matter

LDers0 = OrderedLDersAt_0
Zeta15 = CompNum(exp(2*pi*I/15))
partialZetaVals = []
for i in range(0, 3):
    for j in range(0, 5):
        zetaVal = 0 ##(i,j) - ZetaVal
        for k in range(0, 15):
            zetaVal = zetaVal + zeta3^(-i*k)*zeta5^(-j*k)*LDers0[k]
        zetaVal = (1/15)*zetaVal
        partialZetaVals.append(zetaVal)


## old code order
##for i in range(0, 15):
##    zetaVal = 0 ##(i,i)-ZetaVal
##    for j in range(0, 15):
##        zetaVal = zetaVal + Zeta15^(-i*j)*LDers0[j]
##    zetaVal = (1/15)*zetaVal
##    partialZetaVals.append(zetaVal)

R.<x> = PolynomialRing(QQ)
F.<a> = NumberField(x^2 + 23)
S.<y> = PolynomialRing(F)
H.<b> = F.extension(y^3 - y + 1)
T.<z> = PolynomialRing(H)
H1.<c> = H.extension(z^5 - 10*z^3 + 5*z^2 + 10*z + 1)
U.<w> = PolynomialRing(H1)

## min poly for Stark Units
h = w^15 - 832535*w^14 + 65231675*w^13 - 5650639400*w^12 + 15533478425*w^11 - 39376942640*w^10 - 212804236525*w^9 - 380541320125*w^8 - 2607229594750*w^7 - 2183192838625*w^6 + 3771011381950*w^5 - 1207366794625*w^4 + 99067277500*w^3 - 221569375*w^2 + 466875*w - 125

print "roots"
rootsMults = h.roots()
StarkUnitsExact = []
for i in range(0, 15):
    StarkUnitsExact.append(rootsMults[i][0])

StarkUnitsData = []
for unit in StarkUnitsExact:
    coeffs = []
    for i in range(0, 5):
        for j in range(0, 3):
            for k in range(0, 2):
                coeffs.append(unit.list()[i].list()[j].list()[k])
    StarkUnitsData.append(coeffs)


R.<x> = PolynomialRing(CompNum)
f = x^3 - x + 1

b = CompNum(f.roots()[0][0])
b1 = CompNum(f.roots()[1][0])
b2 = CompNum(f.roots()[2][0])
Hroots = [b, b1, b2]
a = CompNum((-23).sqrt())
c = CompNum(exp(2*pi*I/25) + exp(2*pi*I*7/25) + exp(-2*pi*I/25) + exp(-2*pi*I*7/25))
c1 = CompNum(exp(2*pi*I*2/25) + exp(2*pi*I*11/25) + exp(-2*pi*I*2/25) + exp(-2*pi*I*11/25))
c2 = CompNum(exp(2*pi*I*4/25) + exp(2*pi*I*3/25) + exp(-2*pi*I*4/25) + exp(-2*pi*I*3/25))
c3 = CompNum(exp(2*pi*I*8/25) + exp(2*pi*I*6/25) + exp(-2*pi*I*8/25) + exp(-2*pi*I*6/25))
c4 = CompNum(exp(2*pi*I*9/25) + exp(2*pi*I*12/25) + exp(-2*pi*I*9/25) + exp(-2*pi*I*12/25))
Q1roots = [c, c1, c2, c3, c4]

print "StarkUnitsNum"


StarkUnitsNum = []
for n in range(0, 15):
    StarkUnitNum = CompNum(0)
    data = StarkUnitsData[n]
    for i in range(0, 5):
        for j in range(0, 3):
            for k in range(0, 2):
                StarkUnitNum = StarkUnitNum + CompNum(data[6*i + 2*j + k])*c^i*b^j*a^k
    StarkUnitsNum.append(StarkUnitNum)

StarkLogsNum = []
for i in range(0, 15):
    StarkLogsNum.append((-1)*((StarkUnitsNum[i].abs()).log()))

## find u(0,0)

print "finding Stark Units"
print "should only get one number"

orderedStarkUnitsData = []
orderedStarkUnitsLogs = []
partZetaVal = partialZetaVals[0].real_part()
for i in range(0, 15):
    if ((partZetaVal - StarkLogsNum[i]).abs() < 10^(-100)):
        print i
        u00Data = StarkUnitsData[i]

## embed u00 as other Stark units to get correct order

## the Stark units are ordered as (i,j)

for i in range(0, 3):
    for j in range(0, 5):
        uij = embedStarkUnit(u00Data, a, Hroots[i], Q1roots[j], CompNum)
        for k in range(0, 15):
            if ((uij - StarkUnitsNum[k]).abs() < 10^(-100)):
                print "this should appear exactly once for each i,j"
                print [i,j,k]
                orderedStarkUnitsData.append(StarkUnitsData[k])
                orderedStarkUnitsLogs.append(StarkLogsNum[k])

for i in range(0, 14):
    print ((partialZetaVals[i] - orderedStarkUnitsLogs[i]).abs() < 10^(-100))



























