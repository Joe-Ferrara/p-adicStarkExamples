
C.<zeta25> = CyclotomicField(25)
d = zeta25 - 1
g = (zeta25 - 1).minpoly()

p = 5
prec = 72
ramInd = 20
sizeResField = 25
pAdics = Qp(p, prec)
Rp.<X> = PolynomialRing(pAdics)
g = Rp(g)
CpUnif.<D> = pAdics.extension(g)
Zeta25 = D + 1
Zeta5 = Zeta25^5

## making sqrtNeg3/sqrtNeg23
h = X^2 - pAdics(3/(23))
hroots = [h.roots()[0][0], h.roots()[1][0]]
if (hroots[0].expansion()[0] == 1):
    sqrtNeg3 = hroots[0]
if (hroots[1].expansion()[1] == 1):
    sqrtNeg3 = hroots[1]
print "sqrtNeg3 is"
print sqrtNeg3


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
    ## the coefficients of that element mod p^prec
    data = algNum.list()
    ans = []
    for i in range(0, len(data)):
        temp = data[i] % (5^(72))
        if temp > (5^(72))/2:
            temp = -(5^(72) - temp)
        ans.append(temp)
    return ans



## load and make the LValues, only take the four LValues that we need

LValuesAt_s_0_Data = load("LValuesDataAt_s_0FirstLayer")

print "making L-values"

LValuesAt_s_0 = []
for data in LValuesAt_s_0_Data:
    LValue = make_val_from_data(data)
    LValuesAt_s_0.append(LValue)

## these LValues are ordered for chipsi^i, 1\leq i \leq 4
## must reverse the original order because of the inverse in the integral 

LValues0Pos1 = [LValuesAt_s_0[4], LValuesAt_s_0[3], LValuesAt_s_0[2], LValuesAt_s_0[1]] ## Pos1 refers to alpha = 1 (p-stabilization)

LValues0Pos1Ratios = []
for i in range(0, len(LValues0Pos1)):
    for j in range(i+1, len(LValues0Pos1)):
        ratio = LValues0Pos1[i]/LValues0Pos1[j]
        LValues0Pos1Ratios.append(ratio)


## make the Stark sums and Gauss sums

print "making p-Adic logs"

StLogsPos1Data = load("pAdicLogsPos1Data")

StLogsPos1 = []

for i in range(0, len(StLogsPos1Data)):
    dataPos1 = StLogsPos1Data[i]
    loguPos1 = make_val_from_data(dataPos1)
    StLogsPos1.append(loguPos1)

StSums0Pos1 = [] ## put the Stark sums at 0 here

print "making Stark Sums"

## first we do Pos1

for i in range(1, 5):
    StarkSum = 0
    for j in range(0, 5):
        StarkSum = StarkSum + Zeta5^(i*j)*(StLogsPos1[j] - StLogsPos1[j + 5])
    StSums0Pos1.append(StarkSum)



print "doing ratio stuff"

RatiosPos1_at_0 = [] ## Pos1 corresponds to the p-stabilization

for i in range(0, 4):
    ratio0Pos1 = StSums0Pos1[i]/LValues0Pos1[i] ## since at 0 same alpha on both sides
    RatiosPos1_at_0.append(ratio0Pos1)

GaussSums = []
for k in range(1, 5):
    GaussSum = 0
    for j in range(0, 20):
        GaussSum = (Zeta5^(k*j))*(Zeta25^(2^j)) + GaussSum
    GaussSums.append(GaussSum)
GaussSumRatios = []
for i in range(0, 3):
    for j in range(i + 1, 4):
        GSRatio = GaussSums[3 - j]/GaussSums[3 - i]
        GaussSumRatios.append(GSRatio)

StarkSumRatios0Pos1 = []
for i in range(0, len(StSums0Pos1)):
    for j in range(i + 1, len(StSums0Pos1)):
        ratio = StSums0Pos1[i]/StSums0Pos1[j]
        StarkSumRatios0Pos1.append(ratio)

unitAndGSRatiosPos1 = []
for i in range(0, len(GaussSums)):
    for j in range(i + 1, len(GaussSums)):
        ratio = (StSums0Pos1[i]*GaussSums[3 - i])/(StSums0Pos1[j]*GaussSums[3-j])
        unitAndGSRatiosPos1.append(ratio)
        


RatioCompsPos1_at_0 = []
RatioCompsPos1_at_0_GS = []

for i in range(0, 4):
    for j in range(i + 1, 4):
        RatioCompsPos1_at_0.append(RatiosPos1_at_0[i]/RatiosPos1_at_0[j])
        RatioCompsPos1_at_0_GS.append((RatiosPos1_at_0[i]*GaussSums[3 - i])/(RatiosPos1_at_0[j]*GaussSums[3-j]))

RatioComps = [RatioCompsPos1_at_0]

print "doing approximations"

approximations = []
elementsOfC = []
temp1 = []
temp2 = []
for j in range(0, 6):
    ratio = RatioCompsPos1_at_0[j]
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    temp1.append(alg)
    approx = approximate(alg)
    temp2.append(approx)
elementsOfC.append(temp1)
approximations.append(temp2)  

approximationsGS = []
elementsOfC_GS = []
temp1 = []
temp2 = []
for j in range(0, 6):
    ratio = RatioCompsPos1_at_0_GS[j]
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    temp1.append(alg)
    approx = approximate(alg)
    temp2.append(approx)
elementsOfC_GS.append(temp1)
approximationsGS.append(temp2)  

valuations = []
literalDifference = []
for i in range(0, len(RatioCompsPos1_at_0)):
    diff = RatioCompsPos1_at_0[i] - GaussSumRatios[i]
    literalDifference.append(diff)
    val = (diff).valuation()
    valuations.append(val)

## pi-adic expansion





          


