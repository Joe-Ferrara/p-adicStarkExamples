## check that we get what we predict at s = 1


StLogsPos1Pre = load("pAdicLogsPos1Proj0")
StLogsNeg1Pre = load("pAdicLogsNeg1Proj0")

StLogsPos1Data = []
StLogsNeg1Data = []
for i in range(0, len(StLogsPos1Pre)):
    tempPos1 = StLogsPos1Pre[i]
    tempNeg1 = StLogsNeg1Pre[i]
    dataPos1 = [tempPos1.valuation(), tempPos1.expansion()]
    dataNeg1 = [tempNeg1.valuation(), tempNeg1.expansion()]
    StLogsPos1Data.append(dataPos1)
    StLogsNeg1Data.append(dataNeg1)


C.<zeta25> = CyclotomicField(25)
d = zeta25 - 1
g = (zeta25 - 1).minpoly()

p = 5
prec = 60
ramInd = 20
sizeResField = 25
pAdics = Qp(p, prec)
Rp.<X> = PolynomialRing(pAdics)
g = Rp(g)
CpUnif.<D> = pAdics.extension(g)
Zeta25 = D + 1
Zeta5 = Zeta25^5


## functions

def make_val_from_data(data):
    if len(data[1]) == 0:
        return 0
    else:
        ans = 0
        for i in range(0, len(data[1])):
            ans = CpUnif(data[1][i])*D^(data[0] + i) + ans
        return ans

def make_algebraic(data):
    ## take the p-adic numbers and return elements of C
    ans = 0
    for i in range(0, len(data[1])):
        ans = data[1][i]*d^(data[0] + i) + ans
    return ans

def approximate(algNum):
    ## take an element of C and return
    ## the coefficients of that element mod p^(prec - 10)
    data = algNum.list()
    ans = []
    for i in range(0, len(data)):
        temp = data[i] % (5^(50))
        if temp > (5^(50))/2:
            temp = -(5^(50) - temp)
        ans.append(temp)
    return ans



## load and make the LValues, only take the four LValues that we need

LValuesAt_s_1_Data = load("LValuesDataAt_s_1_new")

LValuesAt_s_1 = []
for data in LValuesAt_s_1_Data:
    LValue = make_val_from_data(data)
    LValuesAt_s_1.append(LValue)

LValuesPos1 = []
for i in range(1, 5):
    LValuesPos1.append(LValuesAt_s_1[i])

LValuesNeg1 = []
for i in range(6, 10):
    LValuesNeg1.append(LValuesAt_s_1[i])


## make the Stark sums and Gauss sums

StLogsPos1 = []
StLogsNeg1 = []

for i in range(0, len(StLogsPos1Data)):
    dataPos1 = StLogsPos1Data[i]
    dataNeg1 = StLogsNeg1Data[i]
    loguPos1 = make_val_from_data(dataPos1)
    StLogsPos1.append(loguPos1)
    loguNeg1 = make_val_from_data(dataNeg1)
    StLogsNeg1.append(loguNeg1)


StSumsPos1 = [] ## put the 4 Stark sums for this embedding here
StSumsNeg1 = []


for k in range(1, 5):
    StSumPos1 = 0
    StSumNeg1 = 0
    for j in range(0, 10):
        StSumPos1 = (CpUnif((-1)^j))*(Zeta5^(k*j))*StLogsPos1[j] + StSumPos1
        StSumNeg1 = (CpUnif((-1)^j))*(Zeta5^(k*j))*StLogsNeg1[j] + StSumNeg1
    StSumsPos1.append(StSumPos1)
    StSumsNeg1.append(StSumNeg1)

## StSumPos1 and StSumNeg1 are the values you want


## make the Gauss sums
GaussSums = []
for k in range(1, 5):
    GaussSum = 0
    for j in range(0, 20):
        GaussSum = (Zeta5^(k*j))*(Zeta25^(2^j)) + GaussSum
    GaussSums.append(GaussSum)

## make ratios
RatiosPos1 = [] ## Pos1 corresponds to the p-stabilization chosen, the StSum will be with Neg1
RatiosNeg1 = [] ## Same convention as above
for k in range(0, 4):
    ratioPos1 = (LValuesPos1[k])/(StSumsNeg1[k])
    ratioNeg1 = (LValuesNeg1[k])/(StSumsPos1[k])
    RatiosPos1.append(ratioPos1)
    RatiosNeg1.append(ratioNeg1)

RatioComparisonsPos1 = []
RatioComparisonsNeg1 = []
for i in range(0, 3):
    for j in range(i + 1, 4):
        RatioPos1 = RatiosPos1[i]/RatiosPos1[j]
        RatioNeg1 = RatiosNeg1[i]/RatiosNeg1[j]
        RatioComparisonsPos1.append(RatioPos1)
        RatioComparisonsNeg1.append(RatioNeg1)

## test algebraicity
approx1Pos1 = []
for ratio in RatioComparisonsPos1:
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Pos1.append(approx)
approx1Neg1 = []
for ratioPair in RatioComparisonsNeg1:
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Neg1.append(approx)

print "approx1Pos1"
print approx1Pos1
print "approx1Neg1"
print approx1Neg1

## check against the Gauss sum ratios
GaussSumRatios = []
for i in range(0, 3):
    for j in range(i + 1, 4):
        GSRatio = GaussSums[j]/GaussSums[i]
        GaussSumRatios.append(GSRatio)

valuationsPos1 = []
valuationsNeg1 = []
for i in range(0, len(RatioComparisonsPos1)):
    val = (RatioComparisonsPos1[i] - GaussSumRatios[i]).valuation()
    valuationsPos1.append(val)
    val = (RatioComparisonsNeg1[i] - GaussSumRatios[i]).valuation()
    valuationsNeg1.append(val)

print "valuationsPos1"
print valuationsPos1
print "valuationsNeg1"
print valuationsNeg1








