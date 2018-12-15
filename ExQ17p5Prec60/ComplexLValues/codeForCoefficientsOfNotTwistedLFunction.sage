load("Functions.sage")

## set parameters

## set parameters

R.<x> = PolynomialRing(QQ)
K.<alpha> = NumberField(x^2 - 17)
S.<y> = PolynomialRing(K)
L.<beta> = K.extension(y^2 - (4 + alpha))

## for first attempt, try 10,000 coefficients

coefficients = []
for n in range(1, 10001):
    coeff = coeffLfunc(n, K, L, 1)
    coefficients.append(coeff)
    save(coefficients, "coeffsNotTwistedLfunc")