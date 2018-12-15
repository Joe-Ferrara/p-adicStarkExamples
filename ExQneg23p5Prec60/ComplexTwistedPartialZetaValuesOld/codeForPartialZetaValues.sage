## Code for making partial zeta values. The code is short, but the ideas
## are somewhat complicated. There should be a paper explanation to go
## with this.

#############################################################################
#############################################################################
##  Note: The element (2, zeta3) in the group                              ##
##               (Z/25Z)^*/<2^5> x {1,zeta3, zeta3^2}                      ##
##        generates the Galois group that we are interested in. The        ##
##        partial zeta values are ordered as (2, zeta3)^j as j runs        ##
##        from 0 to 15. We use the relation between L-functions and        ##
##        partial zeta functions to determine the partial zeta values.     ##
##        We also use the fact that there are 8 characters that cut out    ##
##        the field we are interested in, and only for those 8 characters  ##
##        is the first derivative of the partial zeta functions (and       ##
##        L-functions) not zero. Also, we use that there is a redundancy   ##
##        in the L-functions, which gives pairs of L-functions which are   ##
##        the same (see the paper explanation).                            ##
#############################################################################
#############################################################################

## parameters

prec = 300
C = ComplexField(300)
Cyc.<zeta15> = CyclotomicField(15)
zeta5 = zeta15^3
zeta3 = zeta15^5

## L-values

LDers_at_0 = load("twistedLDers_at_0")

## make the characters

characters = []
G.<psi> = DirichletGroup(25, Cyc, zeta5, 5)
for i in range(1, 5):
    char = psi^i
    characters.append(char)

## make partial zeta values

partialZetaValues = []
for j in range(0, 15):
    partialZetaVal = 0
    for i in range(0, 4):
        charVal = C(((zeta3^2)^j + zeta3^j)*characters[3 - i](2^j))
        LVal = LDers_at_0[i]
        const = 1/15 ## denominator is Galois group order
        partialZetaVal = const*charVal*LVal + partialZetaVal
    partialZetaValues.append(partialZetaVal)

## note that the partial zeta values are real numbers
## the imaginary part is extremely close to zero, so we eliminate it

for i in range(0, 15):
    partialZetaValue = partialZetaValues[i].real_part()
    partialZetaValues[i] = partialZetaValue

## save stuff

save(partialZetaValues, "partialZetaValues")






