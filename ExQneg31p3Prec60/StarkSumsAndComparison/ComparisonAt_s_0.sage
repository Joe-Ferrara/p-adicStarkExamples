

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
        temp = data[i] % (3^(50))
        if temp > (3^(50))/2:
            temp = -(3^(50) - temp)
        ans.append(temp)
    return ans



## load and make the LValues, only take the four LValues that we need

LValuesAt_s_0_Data = load("LValuesDataAt_s_0FirstLayer")

print "making L-values"

LValuesAt_s_0 = []
for data in LValuesAt_s_0_Data:
    LValue = make_val_from_data(data)
    LValuesAt_s_0.append(LValue)

LValuesPos1 = [LValuesAt_s_0[2], LValuesAt_s_0[1]] ## Pos1 refers to alpha = 1 (p-stabilization)

LValuesNeg1 = [LValuesAt_s_0[5], LValuesAt_s_0[4]] ## Neg1 refers to alpha = -1 (p-stabilization)


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


StSumsPos1Try1 = [] ## put the 2 Stark sums here
StSumsNeg1Try1 = []
StSumsPos1Try2 = []
StSumsNeg1Try2 = []

## note the rho-projection not just chi projection

print "making Stark Sums"

for k in range(1, 3):
    StSumPos1Try1 = 0
    StSumNeg1Try1 = 0
    StSumPos1Try2 = 0
    StSumNeg1Try2 = 0
    for j in range(0, 3):
        for i in range(0, 3):
            StSumPos1Try1 = ((Zeta3^i)*(Zeta3^(k*j)) + (Zeta3^(2*i))*(Zeta3^(k*j)))*StLogsPos1[3*j + i] + StSumPos1Try1
            StSumPos1Try2 = ((Zeta3^i)*(Zeta3^(k*j)))*StLogsPos1[3*j + i] + StSumPos1Try2
            StSumNeg1Try1 = ((Zeta3^i)*(Zeta3^(k*j)) + (Zeta3^(2*i))*(Zeta3^(k*j)))*StLogsNeg1[3*j + i] + StSumNeg1Try1
            StSumNeg1Try2 = ((Zeta3^i)*(Zeta3^(k*j)))*StLogsNeg1[3*j + i] + StSumNeg1Try2
    StSumsPos1Try1.append(StSumPos1Try1)
    StSumsNeg1Try1.append(StSumNeg1Try1)
    StSumsPos1Try2.append(StSumPos1Try2)
    StSumsNeg1Try2.append(StSumNeg1Try2)


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

RatiosPos1Try1 = [] ## Pos1 corresponds to the p-stabilization chosen, the StSum will be with Pos1
##RatiosNeg1Try1 = [] ## Same convention as above
############################################################
## THE RATIO MAY OR MAY NOT BE RIGHT - TRIED BOTH OF THEM ##
############################################################
for k in range(0, 2):
    ratioPos1 = (LValuesPos1[k])/(StSumsPos1Try1[k])
##    ratioNeg1 = (LValuesNeg1[k])/(StSumsNeg1Try1[k])
    RatiosPos1Try1.append(ratioPos1)
##    RatiosNeg1Try1.append(ratioNeg1)


RatiosPos1Try2 = [] ## Pos1 corresponds to the p-stabilization chosen, the StSum will be with Pos1
RatiosNeg1Try2 = [] ## Same convention as above
## Note that because at s=0 the order is reversed compared to s=1
for k in range(0, 2):
    ratioPos1 = (LValuesPos1[k])/(StSumsPos1Try2[k])
    ratioNeg1 = (LValuesNeg1[k])/(StSumsNeg1Try2[k])
    RatiosPos1Try2.append(ratioPos1)
    RatiosNeg1Try2.append(ratioNeg1)



RatioComparisonsPos1Try1 = [RatiosPos1Try1[0]/RatiosPos1Try1[1]]
##RatioComparisonsNeg1Try1 = [RatiosNeg1Try1[0]/RatiosNeg1Try1[1]]


RatioComparisonsPos1Try2 = [RatiosPos1Try2[0]/RatiosPos1Try2[1]]
RatioComparisonsNeg1Try2 = [RatiosNeg1Try2[0]/RatiosNeg1Try2[1]]


## check algebraicity
approx1Pos1Try1 = []
for ratio in RatioComparisonsPos1Try1:
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Pos1Try1.append(approx)
##approx1Neg1Try1 = []
##for ratio in RatioComparisonsNeg1Try1:
##    data = [ratio.valuation(), ratio.expansion()]
##    alg = make_algebraic(data)
##    approx = approximate(alg)
##    approx1Neg1Try1.append(approx)

## check algebraicity
approx1Pos1Try2 = []
for ratio in RatioComparisonsPos1Try2:
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Pos1Try2.append(approx)
approx1Neg1Try2 = []
for ratio in RatioComparisonsNeg1Try2:
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Neg1Try2.append(approx)



## compare to the GaussSumRatios
GaussSumRatios = [GaussSums[1]/GaussSums[0]]

valuationsPos1Try1 = []
##valuationsNeg1Try1 = []
for i in range(0, len(RatioComparisonsPos1Try1)):
    val = (RatioComparisonsPos1Try1[i] - GaussSumRatios[i]).valuation()
    valuationsPos1Try1.append(val)
##    val = (RatioComparisonsNeg1Try1[i] - GaussSumRatios[i]).valuation()
##    valuationsNeg1Try1.append(val)

print "valuationsPos1Try1"
print  valuationsPos1Try1
##print "valuationsNeg1Try1"
##print  valuationsNeg1Try1
    
valuationsPos1Try2 = []
valuationsNeg1Try2 = []
for i in range(0, len(RatioComparisonsPos1Try2)):
    val = (RatioComparisonsPos1Try2[i] - GaussSumRatios[i]).valuation()
    valuationsPos1Try2.append(val)
    val = (RatioComparisonsNeg1Try2[i] - GaussSumRatios[i]).valuation()
    valuationsNeg1Try2.append(val)

print "valuationsPos1Try2"
print  valuationsPos1Try2
print "valuationsNeg1Try2"
print  valuationsNeg1Try2
    


