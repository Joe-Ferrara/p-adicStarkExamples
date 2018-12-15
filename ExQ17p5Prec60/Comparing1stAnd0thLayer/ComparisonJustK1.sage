print 1

StLogsPos1Pre = load("pAdicLogsPos1Proj0")
StLogsNeg1Pre = load("pAdicLogsNeg1Proj0")

print 2

StLogsPos1Data = []
StLogsNeg1Data = []
for i in range(0, len(StLogsPos1Pre)):
    tempPos1 = StLogsPos1Pre[i]
    tempNeg1 = StLogsNeg1Pre[i]
    dataPos1 = [tempPos1.valuation(), tempPos1.expansion()]
    dataNeg1 = [tempNeg1.valuation(), tempNeg1.expansion()]
    StLogsPos1Data.append(dataPos1)
    StLogsNeg1Data.append(dataNeg1)

print 3

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

print 4

## functions

def make_val_from_data(data):
    if len(data[1]) == 0:
        return 0
    else:
        ans = 0
        for i in range(0, len(data[1])):
            ans = CpUnif(data[1][i])*D^(data[0] + i) + ans
        return ans

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

print 5

for k in range(1, 5):
    StSumPos1 = 0
    StSumNeg1 = 0
    for j in range(0, 10):
        ## these are the Stark sums form chi*psi^k at s = 1 (notice the -k in the exponent)
        StSumPos1 = (CpUnif((-1)^j))*((Zeta5^(-k))^j)*StLogsPos1[j] + StSumPos1
        StSumNeg1 = (CpUnif((-1)^j))*((Zeta5^(-k))^j)*StLogsNeg1[j] + StSumNeg1
    StSumsPos1.append(StSumPos1)
    StSumsNeg1.append(StSumNeg1)

## StSumPos1 and StSumNeg1 are the values you want

print 6

## make the Gauss sum
GaussSums = []
for k in range(1, 5):
    GaussSum = 0
    for j in range(0, 20):
        ## this is the Gauss sum for psi^k
        GaussSum = ((Zeta5^k)^j)*(Zeta25^(2^j)) + GaussSum
    GaussSums.append(GaussSum)

## make tau(psi^k)*StSum(psi^k) which should be in Q_p(zeta5)

print 7

GSStSumsPos1 = []
GSStSumsNeg1 = []
for k in range(0, 4):
    GSStSumPos1 = GaussSums[k]*StSumsPos1[k]
    GSStSumsPos1.append(GSStSumPos1)
    GSStSumNeg1 = GaussSums[k]*StSumsNeg1[k]
    GSStSumsNeg1.append(GSStSumNeg1)


## now make the LValues

LValuesAt_s_1_Data = load("LValuesDataAt_s_1_new")

## the LValues are ordered as Lp(phi_alpha, psi^(-k), 1) for alpha = 1, k = 0,1,2,3,4 then alpha = -1, k = 0,1,2,3,4
## here Lp(phi_alpha, psi^(-k), 1) = int_(Z_p^*) psi^k dmu_alpha where mu_alpha = phi_alpha(0 - infty)

print 8

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

## compare the ratios Lp(phi_alpha, psi^i, 1)/Lp(phi_alpha, psi^j, 1) and tau(psi^i)logp(u_chi*psi^i)/tau(psi^j)logp(u_chi*psi^j)

print 9

LValRatiosPos1 = []
StSumRatiosPos1 = []
LValRatiosNeg1 = []
StSumRatiosNeg1 = []
for i in range(0, 4):
    for j in range(i + 1, 4):
        LValRatio = LValuesPos1[3 - i]/LValuesPos1[3 - j]
        LValRatiosPos1.append(LValRatio)
        LValRatio = LValuesNeg1[3 - i]/LValuesNeg1[3 - j]
        LValRatiosNeg1.append(LValRatio)
        StSumRatio = GSStSumsPos1[i]/GSStSumsPos1[j]
        StSumRatiosPos1.append(StSumRatio)
        StSumRatio = GSStSumsNeg1[i]/GSStSumsNeg1[j]
        StSumRatiosNeg1.append(StSumRatio)

print 10

valuationsPos1 = []
valuationsNeg1 = []
for k in range(0, 6):
    val = (LValRatiosNeg1[k] - StSumRatiosPos1[k]).valuation()
    valuationsPos1.append(val)
    val = (LValRatiosPos1[k] - StSumRatiosNeg1[k]).valuation()
    valuationsNeg1.append(val)


## compare the ratios Lp(phi_alpha, psi^i, 1)/Lp(phi_alpha, psi^j, 1) and tau(psi^i)logp(u_chi*psi^i)/tau(psi^j)logp(u_chi*psi^j)

##print 11

##LValRatiosTPos1 = []
##LValRatiosTNeg1 = []
##for i in range(0, 4):
##    for j in range(i + 1, 4):
##        LValRatio = LValuesPos1[3 - i]/LValuesPos1[3 - j]
##        LValRatiosTPos1.append(LValRatio)
##        LValRatio = LValuesNeg1[3 - i]/LValuesNeg1[3 - j]
##        LValRatiosTNeg1.append(LValRatio)

##print 12

##valuationsTPos1 = []
##valuationsTNeg1 = []
##for k in range(0, 6):
##    val = (LValRatiosTNeg1[k] - StSumRatiosPos1[k]).valuation()
##    valuationsTPos1.append(val)
##    val = (LValRatiosTPos1[k] - StSumRatiosNeg1[k]).valuation()
##    valuationsTNeg1.append(val)

print "valuationsPos1"
print valuationsPos1

print "valuationsNeg1"
print valuationsNeg1

##print "valuationsTPos1"
##print valuationsTPos1

##print "valuationsTNeg1"
##print valuationsTNeg1





