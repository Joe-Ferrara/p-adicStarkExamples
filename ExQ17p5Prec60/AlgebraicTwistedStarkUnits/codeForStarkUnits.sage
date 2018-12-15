## This follows the method outlined by Stark in conference proceedings
## Almost an identical method is used and explained in Dummit, Sands, and Tangedal
## See the paper explanation for what is going on and notation

## load

partialZetaValues = load("partialZetaValues")

## functions

def make_bounds(thetaNum, n, d, Reals):
    if d%4 == 1:
        lowerBound = Reals((thetaNum - n)/(d.sqrt())).floor()
        upperBound = Reals((thetaNum + n)/(d.sqrt())).ceiling() + 1
        return [lowerBound, upperBound]
    else:
        lowerBound = Reals((thetaNum - n)/(2*(d.sqrt()))).floor()
        upperBound = Reals((thetaNum + n)/(2*(d.sqrt()))).ceiling() + 1
        return [lowerBound, upperBound]

def make_aNum(thetaNum, d, b, Reals):
    if d%4 == 1:
        aNum = Reals(2*thetaNum - b*(d.sqrt()))
        return aNum
    else:
        aNum = Reals(thetaNum - b*(d.sqrt()))
        return aNum

def alg_coeff_theta(thetaNum, j, sizeOfH, d, Reals):
    n = binomial(sizeOfH, j)*2^j
    bounds = make_bounds(thetaNum, n, d, Reals)
    lowerBound = bounds[0]
    upperBound = bounds[1]
    for b in range(lowerBound, upperBound):
        aNum = make_aNum(thetaNum, d, b, Reals)
        if Reals((aNum.round() - aNum).abs()) < Reals(0.00000000001):
            a = aNum.round()
            return [a,b]

def make_St_units_from_gRoot(eta):
    u = (eta + (eta^2 - 4).sqrt())/2
    ubar = (eta - (eta^2 - 4).sqrt())/2
    return [u, ubar]

## parameters

C = ComplexField(300)
R = RealField(300)
T.<X> = PolynomialRing(R)

S.<x> = PolynomialRing(QQ)
K.<alpha> = NumberField(x^2 - 17)

Cyc.<zeta25> = CyclotomicField(25)
g = (zeta25 - 1).minpoly()
CycUnif.<dekka> = NumberField(g)
U.<y> = PolynomialRing(CycUnif)
M.<beta> = CycUnif.extension(y^4 - 8*y^2 - 1)
alpha = beta^2 - 4
V.<z> = PolynomialRing(M)

## make g numerically

unitsAbs = [C(exp(-2*x)) for x in partialZetaValues]
    ## absolute values of the Stark units
    ## we may assume all the Stark units are positive
    ## therefor unitsAbs in the Stark units numerically

gRootsNum = []
for i in range(0, 5):
    eta = unitsAbs[i] + unitsAbs[i + 5]
    gRootsNum.append(eta)
linFactors = [X - eta for eta in gRootsNum]
gNum = prod(linFactors)
thetas = gNum.coefficients()

## the thetas are elements of K, so we use assumptions on the
## Stark units to figure out exactly what the thetas are
## that gives us g, then factoring g gives us the etas
## easy to get the Stark units from the etas
## note: here we are using etas where the written
##       explanation has alphas

thetasExact = []
for j in range(0, len(thetas)):
    sizeOfH = 5 ## SPECIFIC FOR THIS EXAMPLE
    thetaNum = thetas[sizeOfH - j]
    d = 17 ## SPECIFIC FOR THIS EXAMPLE
    algCoeffs = alg_coeff_theta(thetaNum, j, sizeOfH, d, R)
    thetaAlg = (algCoeffs[0] + algCoeffs[1]*alpha)/2
    ## NOTE: I AM USING THE d = 17 IMPLICITELY HERE IN ABOVE LINE
    thetasExact.append(thetaAlg)

gEx = 0
for i in range(0, len(thetasExact)):
    sizeOfH = 5 ## SPECIFIC FOR THIS EXAMPLE
    gEx = thetasExact[i]*z^(sizeOfH - i) + gEx
gRoots = gEx.roots()

StarkUnits = []
for root in gRoots:
    eta = root[0]
    units = make_St_units_from_gRoot(eta)
    StarkUnits.append(units[0])
    StarkUnits.append(units[1])

## save stuff

## save(StarkUnits, 'algebraicStarkUnits')
