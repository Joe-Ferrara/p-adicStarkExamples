p = 5
prec = 60
Q5 = Qp(p, prec)
R.<X> = PolynomialRing(Q5)
F5.<a> = Q5.extension(X^2 - 17)

f = X^2 + 1
sqrtNeg1 = f.roots()[0][0]
## pick the square root of -1 that is 2 mod 5
if Integers(5)(sqrtNeg1) == Integers(5)(3):
    sqrtNeg1 = f.roots()[1][0]
f = X^2 - (2 - (1/2)*sqrtNeg1)
const0 = f.roots()[0][0]
## pick the root that is 1 mod 5
if Integers(5)(const0) == Integers(5)(4):
    const0 = f.roots()[1][0]
const1 = 1/(2*const0)
## B is sqrt(4 + sqrt(17)) a root of x^2 - (4 + sqrt(17))
B = const0 + const1*a
Bbar = const0 - const1*a

## this is Stark unit in K coming from the TwoDifferentRepresentation... sage worksheet
StarkUnit = (((1/4)*a - 3/4)*B + (1/4)*a + 1/4)^(-4)

## get conjugate of Stark unit
StarkUnitBar = ((-(1/4)*a - 3/4)*Bbar - (1/4)*a + 1/4)^(-4)

## get the projections
StarkUnitPos1 = StarkUnit*StarkUnitBar
StarkUnitNeg1 = StarkUnit/StarkUnitBar

ploguPos1 = StarkUnitPos1.log()
ploguPos1Data = [ploguPos1.valuation(), ploguPos1.list()]
save(ploguPos1Data, 'ploguPos1KunitData')
ploguNeg1 = StarkUnitNeg1.log()/a
ploguNeg1Data = [ploguNeg1.valuation(), ploguNeg1.list()]
save(ploguNeg1Data, 'ploguNeg1KunitData')

print ploguPos1
print ploguNeg1
