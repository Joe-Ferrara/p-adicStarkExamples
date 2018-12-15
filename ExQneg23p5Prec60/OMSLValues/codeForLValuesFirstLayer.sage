## note that we've implicitely chosen am embedding here for zeta25, the 25th primitive root of unity

## make sure you've loaded master.sage

eigenVals = [1,-1]

print "loading plus eigensymbols"
eigenSymbs = load("plusHeckeSymbs")


print "making the LValues at s = 1"
LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = getCosetMeasures(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 5):
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_1FirstLayer')

print "making the data for LValues at s = 1"
LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.expansion()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_1FirstLayer')

## do the LValue at 0

eigenVals = [1, -1]

print "loading the minus eigensymbolss"
eigenSymbs = load('minusHeckeSymbs')

print "making the LValues at s = 0"
LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = cosetIntegralsAt_s_0(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 5):
        print j
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_0FirstLayer')

print "making the data for the LValues"
LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.expansion()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_0FirstLayer')
