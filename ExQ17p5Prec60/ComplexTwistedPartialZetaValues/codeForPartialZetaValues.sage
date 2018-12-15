## Code for making partial zeta values. The code is short, but the ideas
## are somewhat complicated. There should be a paper explanation to go
## with this.

#############################################################################
#############################################################################
##  Note: The element (2, -1) in the group (Z/25Z)^*/<2^5> x {1,-1}        ##
##        generates the Galois group that we are interested in. The        ##
##        partial zeta values are ordered as (2, -1)^j as j runs from 0    ##
##        to 9. We use the relation between L-functions and partial zeta   ##
##        functions to determine the partial zeta values. We also use the  ##
##        fact that there are four characters, that cut out the field we   ##
##        are interested in, and only for those four characters, is the    ##
##        first derivative of the partial zeta functions (and L-functions) ##
##        not zero.                                                        ##
#############################################################################
#############################################################################

## parameters

prec = 300
C = ComplexField(300)
Cyc.<zeta5> = CyclotomicField(5)

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
for j in range(0, 10):
    partialZetaVal = 0
    for i in range(0, 4):
        charVal = C((-1)^j*characters[3 - i](2^j))
        LVal = LDers_at_0[i]
        const = 1/10 ## denominator is Galois group order
        partialZetaVal = const*charVal*LVal + partialZetaVal
    partialZetaValues.append(partialZetaVal)

## note that the partial zeta values are real numbers
## the imaginary part is extremely close to zero, so we eliminate it

for i in range(0, 10):
    partialZetaValue = partialZetaValues[i].real_part()
    partialZetaValues[i] = partialZetaValue

## save stuff

save(partialZetaValues, "partialZetaValues")






