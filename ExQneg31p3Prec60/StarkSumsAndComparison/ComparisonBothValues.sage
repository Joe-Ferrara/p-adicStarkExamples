

C.<zeta9> = CyclotomicField(9)
d = zeta9 - 1
g = (zeta9 - 1).minpoly()

p = 3
prec = 60
ramInd = 6
sizeResField = 9
pAdics = Qp(p, prec)
Rp.<X> = PolynomialRing(pAdics)
g = Rp(g)
CpUnif.<D> = pAdics.extension(g)
Zeta9 = D + 1
Zeta3 = Zeta9^3


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
        temp = data[i] % (3^(40))
        if temp > (3^(40))/2:
            temp = -(3^(40) - temp)
        ans.append(temp)
    return ans



## load and make the LValues, only take the four LValues that we need

LValuesAt_s_0_Data = load("LValuesDataAt_s_0FirstLayer")
LValuesAt_s_1_Data = load("LValuesDataAt_s_1FirstLayer")

print "making L-values"

LValuesAt_s_0 = []
for data in LValuesAt_s_0_Data:
    LValue = make_val_from_data(data)
    LValuesAt_s_0.append(LValue)

LValues0Pos1 = [LValuesAt_s_0[2], LValuesAt_s_0[1]] ## Pos1 refers to alpha = 1 (p-stabilization)

LValues0Neg1 = [LValuesAt_s_0[5], LValuesAt_s_0[4]] ## Neg1 refers to alpha = -1 (p-stabilization)

print "making L-values"

LValuesAt_s_1 = []
for data in LValuesAt_s_1_Data:
    LValue = make_val_from_data(data)
    LValuesAt_s_1.append(LValue)

LValues1Pos1 = [LValuesAt_s_1[2], LValuesAt_s_1[1]] ## Pos1 refers to alpha = 1 (p-stabilization)

LValues1Neg1 = [LValuesAt_s_1[5], LValuesAt_s_1[4]] ## Neg1 refers to alpha = -1 (p-stabilization)


## make the Stark sums and Gauss sums

print "making p-Adic logs"

StLogsPos1Data = load("pAdicLogsPos1Data")
StLogsNeg1Data = load("pAdicLogsNeg1Data")

StLogsPos1 = []
StLogsNeg1 = []

for i in range(0, len(StLogsPos1Data)):
    dataPos1 = StLogsPos1Data[i]
    dataNeg1 = StLogsNeg1Data[i]
    loguPos1 = make_val_from_data(dataPos1)
    StLogsPos1.append(loguPos1)
    loguNeg1 = make_val_from_data(dataNeg1)
    StLogsNeg1.append(loguNeg1)


StSumsPos1 = [] ## put the 2 Stark sums here
StSumsNeg1 = []

## doing the rho, chi, and chi^2 projections

print "making Stark Sums"

for k in range(1, 3):
    StSumPos1rho = 0
    StSumPos1chi = 0
    StSumPos1chi2 = 0
    StSumNeg1rho = 0
    StSumNeg1chi = 0
    StSumNeg1chi2 = 0
    for j in range(0, 3):
        for i in range(0, 3):
            StSumPos1rho = (Zeta3^i + Zeta3^(2*i))*(Zeta3^(k*j))*StLogsPos1[3*j + i] + StSumPos1rho
            StSumNeg1rho = (Zeta3^i + Zeta3^(2*i))*(Zeta3^(k*j))*StLogsNeg1[3*j + i] + StSumNeg1rho ## should be 0
            StSumPos1chi = (Zeta3^i)*(Zeta3^(k*j))*StLogsPos1[3*j + i] + StSumPos1chi
            StSumNeg1chi = (Zeta3^i)*(Zeta3^(k*j))*StLogsNeg1[3*j + i] + StSumNeg1chi
            StSumPos1chi2 = (Zeta3^(2*i))*(Zeta3^(k*j))*StLogsPos1[3*j + i] + StSumPos1chi2 ## should be same as StSumPos1chi
            StSumNeg1chi2 = (Zeta3^(2*i))*(Zeta3^(k*j))*StLogsNeg1[3*j + i] + StSumNeg1chi2 ## should be the negative of StSumNeg1chi
    StSumsPos1.append([StSumPos1rho, StSumPos1chi, StSumPos1chi2])
    StSumsNeg1.append([StSumNeg1rho, StSumNeg1chi, StSumNeg1chi2])


## make the Gauss sum

print "making Gauss sums"

GaussSums = []
for k in range(1, 3):
    GaussSum = 0
    for j in range(0, 6):
        GaussSum = (Zeta3^(k*j))*(Zeta9^(2^j)) + GaussSum
    GaussSums.append(GaussSum)

print "doing ratio stuff"

## try 1

RatiosPos1_at_0 = [] ## Pos1 corresponds to the p-stabilization chosen
RatiosNeg1_at_0 = [] ## Same convention as above
for k in range(0, 2):
    ratioPos1rho = (LValues0Pos1[k])/(StSumsPos1[k][0])
    ratioPos1chi = (LValues0Pos1[k])/(StSumsPos1[k][1])
    ratioPos1chi2 = (LValues0Pos1[k])/(StSumsPos1[k][2])
    RatiosPos1_at_0.append([ratioPos1rho, ratioPos1chi, ratioPos1chi2])
    ## ratioNeg1rho = (LValues0Neg1[k])/(StSumsNeg1[k][0]) THE STARK SUM IS 0
    ratioNeg1chi = (LValues0Neg1[k])/(StSumsNeg1[k][1])
    ratioNeg1chi2 = (LValues0Neg1[k])/(StSumsNeg1[k][2])
    RatiosNeg1_at_0.append([CpUnif(1), ratioNeg1chi, ratioNeg1chi2])

RatiosPos1_at_1 = []
RatiosNeg1_at_1 = []
for k in range(0, 2):
    ## ratioPos1rho = (LValues1Pos1[k])/(StSumsNeg1[1-k][0]) THE STARK SUM IS 0
    ratioPos1chi = (LValues1Pos1[k])/(StSumsNeg1[1-k][2])
    ratioPos1chi2 = (LValues1Pos1[k])/(StSumsNeg1[1-k][1])
    RatiosPos1_at_1.append([CpUnif(1), ratioPos1chi, ratioPos1chi2])
    ratioNeg1rho = (LValues1Neg1[k])/(StSumsPos1[1-k][0])
    ratioNeg1chi = (LValues1Neg1[k])/(StSumsPos1[1-k][2])
    ratioNeg1chi2 = (LValues1Neg1[k])/(StSumsPos1[1-k][1])
    RatiosNeg1_at_1.append([ratioNeg1rho, ratioNeg1chi, ratioNeg1chi2])


RatioCompsPos1_at_0 = []
RatioCompsNeg1_at_0 = []
RatioCompsPos1_at_1 = []
RatioCompsNeg1_at_1 = []
for i in range(0, 3):
    RatioCompsPos1_at_0.append(RatiosPos1_at_0[0][i]/RatiosPos1_at_0[1][i])
    RatioCompsNeg1_at_0.append(RatiosNeg1_at_0[0][i]/RatiosNeg1_at_0[1][i])
    RatioCompsPos1_at_1.append(RatiosPos1_at_1[0][i]/RatiosPos1_at_1[1][i])
    RatioCompsNeg1_at_1.append(RatiosNeg1_at_1[0][i]/RatiosNeg1_at_1[1][i])

RatioComps = [RatioCompsPos1_at_0, RatioCompsNeg1_at_0, RatioCompsPos1_at_1, RatioCompsNeg1_at_1]

approximations = []
elementsOfC = []
for i in range(0, 4):
    temp1 = []
    temp2 = []
    for j in range(0, 3):
        ratio = RatioComps[i][j]
        data = [ratio.valuation(), ratio.expansion()]
        alg = make_algebraic(data)
        temp1.append(alg)
        approx = approximate(alg)
        temp2.append(approx)
    elementsOfC.append(temp1)
    approximations.append(temp2)   

## Divide by Gauss Sums

GaussSumRatio = GaussSums[1]/GaussSums[0]

RatioCompsGS = [RatioCompsPos1_at_0[1]/GaussSumRatio, RatioCompsNeg1_at_0[1]/GaussSumRatio]
approximationsWGS = []
for i in range(0, 2):
    ratio = RatioCompsGS[i]
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approximationsWGS.append(approx)


