from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding

## parameters

C = ComplexField(300)

## functions

def embedding(element, ComplexField):
    sqrt17 = ComplexField(17.sqrt())
    betaNum = ComplexField((4 + sqrt17).sqrt())
    dekkaNum = ComplexField(exp(2*pi*I/25) - 1)
    ans = 0
    for i in range(0, 4):
        for j in range(0, 20):
            coeffCycUnif = element.list()[i]
            coeffQ = coeffCycUnif.list()[j]
            coeffQ = ComplexField(coeffQ)
            ans = coeffQ*betaNum^i*dekkaNum^j + ans
    return ans

## load

algStarkUnits = load("algebraicStarkUnits")
partialZetaValues = load("partialZetaValues")

## check and order Stark units

complexLogStarkUnits = []
for u in algStarkUnits:
    print "embeddingUnit"
    uNum = embedding(u, C)
    loguNum = C((-1/2)*log(uNum))
    complexLogStarkUnits.append(loguNum)

orderedStarkUnits = []
order = []
for zetaVal in partialZetaValues:
    for i in range(0, len(complexLogStarkUnits)):
        print i, "comparison"
        logu = complexLogStarkUnits[i]
        if (zetaVal - logu).abs() < C(0.00000000001):
            orderedStarkUnits.append(algStarkUnits[i])
            order.append(i)

## save stuff

save(orderedStarkUnits, 'orderedStarkUnits')
save(order, 'order')





