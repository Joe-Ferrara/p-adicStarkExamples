## if you have the eigensymbols, this code will make the LValues

## note must change the 5's to 3's!!!!

## make sure you've loaded master.sage from the other directory

## note that I'm going up to p^2 conductor here so to the first layer
## for the initial calulcations only use the first layer
## using only the first layer just gives one ratio
## in the next code when picking out the LValues to use must be careful

eigenVals = [1, -1]

print "loading plus eigensymbols"
eigenSymbs = load("heckePlusEigenSymbs.sobj") ## load the plus eigensymbols

print "making the LVAlues at s = 1"
LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = getCosetMeasures(eigenSymbs[i], eigenVals[i], 2)
    for j in range(0, 3):
        LValue = getLValue(eigenSymbs[i], cosetIntegrals, 2, j)
        LValues.append(LValue)
save(LValues, 'LValuesAt_s_1FirstLayer')

print "making the data for the LValues"
LValuesData = []
for LVal in LValues:
    val =  LVal.valuation()
    coeffs = LVal.expansion()
    LValuesData.append([val, coeffs])
save(LValuesData, 'LValuesDataAt_s_1FirstLayer')

print "loading the minus eigensymbols"
eigenSymbs = load("heckeMinusEigenSymbs.sobj") ## load the minus eigensymbols

print "making the LValues at s = 0"
LValues = []
for i in range(0, len(eigenSymbs)):
    cosetIntegrals = cosetIntegralsAt_s_0(eigenSymbs[i], eigenVals[i], 2)
    print i
    for j in range(0, 3):
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


