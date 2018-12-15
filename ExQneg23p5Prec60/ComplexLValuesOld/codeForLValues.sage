load("Functions.sage")

## some necessary function to clean the code up

def make_coeffs_pairs(w):
    ## hardwired thing
    ## pairs the coefficients of an L-function
    ## with the coefficients of its dual L-function
    ## needs to be done for Dokchitser's code
    ans = [[w[0], w[0]], [w[1], w[4]], [w[2], w[3]], [w[3], w[2]], [w[4], w[1]]]
    return ans

def make_characters(G):
    ans = []
    for psi in G:
        ans.append(psi)
    psi = G.gen()
    char1 = psi^5
    H = DirichletGroup(1)
    char1 = H(char1)
    ans[0] = char1
    return ans

## set parameters

R.<x> = PolynomialRing(QQ)
K.<alpha> = NumberField(x^2 + 23)
S.<y> = PolynomialRing(K)
L.<beta> = K.extension(y^3 - y^2 + 1)
C.<zeta15> = CyclotomicField(15)
zeta5 = zeta15^3
G = DirichletGroup(25, C, zeta5, 5)
characters = make_characters(G)
ComplexNumbers = ComplexField(500)
roots = load("rootNumbersModified.sobj")
rootsComplex = [ComplexNumbers(x) for x in roots]

## make twisted coefficients
## NOTE: the not twisted coeeficients are included
##       as the first element of this list

coefficients = load("coeffsNotTwistedLfunc")
twistedCoeffs = []
for psi in characters:
    psiCoeffs = []
    for n in range(1, 10001):
        coeff = coefficients[n-1]
        newCoeff = C(psi(n))*C(coeff)
        psiCoeffs.append(newCoeff)
    twistedCoeffs.append(psiCoeffs)
save(twistedCoeffs, "twistedCoeffsExact")

## make coefficients complex numbers

for i in range(0, 5):
    coeffs = twistedCoeffs[i]
    coeffs = [ComplexNumbers(coef) for coef in coeffs]
    twistedCoeffs[i] = coeffs
save(twistedCoeffs, "twistedCoeffsComplex")

## pair the coefficients with the inverse character coefficients

twistedCoeffsPairs = make_coeffs_pairs(twistedCoeffs)

## make the conductors

conductors = [23]
for i in range(0, 4):
    conductors.append(23*25*25)

## make Lfunctions using Dokchitser

Lfunctions = []
LfuncCheckFunctEq = []
LValues_at_1 = []
for i in range(0, 5):

    ## set parameters for Dokchitser

    prec = 300
    cond = conductors[i]
    gammaV = [0, 1]
    wt = 1
    poles = []
    residues = 'automatic'
    init = '1'
    coeffs = twistedCoeffsPairs[i]
    root = rootsComplex[i]

    ## make the Lfunctions

    Lfunct = Dokchitser(cond, gammaV, wt, root, poles, residues, prec, init)
    Lfunct.init_coeffs(coeffs[0], 1, coeffs[1])
    Lfunctions.append(Lfunct)
    LfuncCheckFunctEq.append(Lfunct.check_functional_equation())
    LValues_at_1.append(Lfunct(1))

## save stuff

save(Lfunctions, "Lfunctions")
save(LfuncCheckFunctEq, "LfuncCheckFunctEq")
save(LValues_at_1, "LValues_at_1")

## make the derivatives at 0

ComplexNumbers = ComplexField(prec)
LDerivatives_at_0 = []
for i in range(0, 5):
    constant = ComplexNumbers(rootsComplex[i]*(conductors[i].sqrt())*(1/(2*pi)))
    if i == 0:
        LValue_at_1 = LValues_at_1[0]
    else:
        LValue_at_1 = LValues_at_1[5 - i]
    LDerivative_at_0 = constant*LValue_at_1
    LDerivatives_at_0.append(LDerivative_at_0)
save(LDerivatives_at_0, "LDerivatives_at_0") ## save stuff


