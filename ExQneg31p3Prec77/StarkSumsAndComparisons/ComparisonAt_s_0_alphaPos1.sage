

C.<zeta9> = CyclotomicField(9)
d = zeta9 - 1
g = (zeta9 - 1).minpoly()

p = 3
prec = 77
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
    ## the coefficients of that element mod p^prec
    data = algNum.list()
    ans = []
    for i in range(0, len(data)):
        temp = data[i] % (3^(77))
        if temp > (3^(77))/2:
            temp = -(3^(77) - temp)
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

LValuesRatioPos1 = LValuesAt_s_0[2]/LValuesAt_s_0[1]


## make the Stark sums and Gauss sums

print "making p-Adic logs"

StLogsPos1Data = load("pAdicLogsPos1Data")


StLogsPos1 = []

for i in range(0, len(StLogsPos1Data)):
    dataPos1 = StLogsPos1Data[i]
    loguPos1 = make_val_from_data(dataPos1)
    StLogsPos1.append(loguPos1)


StSumsPos1 = []

## note the rho-projection not just chi projection

print "making Stark Sums"

for k in range(1, 3):
    StSumPos1 = 0
    for j in range(0, 3):
        for i in range(0, 3):
            StSumPos1 = ((Zeta3^i)*(Zeta3^(k*j)))*StLogsPos1[3*j + i] + StSumPos1
    StSumsPos1.append(StSumPos1)

## make the Gauss sum

print "making Gauss sums"

GaussSums = []
for k in range(1, 3):
    GaussSum = 0
    for j in range(0, 6):
        GaussSum = (Zeta3^(k*j))*(Zeta9^(2^j)) + GaussSum
    GaussSums.append(GaussSum)

print "doing ratio stuff"


RatiosPos1 = [] ## Pos1 corresponds to the p-stabilization chosen, the StSum will be with Pos1
## Note that because at s=0 the order is reversed compared to s=1
for k in range(0, 2):
    ratioPos1 = (LValuesPos1[k])/(StSumsPos1[k])
    RatiosPos1.append(ratioPos1)

unitAndGSRatioPos1 = (StSumsPos1[0]*GaussSums[1])/(StSumsPos1[1]*GaussSums[0])

RatioComparisonsPos1 = [RatiosPos1[0]/RatiosPos1[1]]

RatioComparisonsPos1WithGS = [(RatiosPos1[0]*GaussSums[0])/(RatiosPos1[1]*GaussSums[1])]

## check algebraicity
algQuantitiesPos1 = []
approx1Pos1 = []
for ratio in RatioComparisonsPos1:
    data = [ratio.valuation(), ratio.expansion()]
    alg = make_algebraic(data)
    algQuantitiesPos1.append(alg)
    approx = approximate(alg)
    approx1Pos1.append(approx)

## compare to the GaussSumRatios
GaussSumRatios = [GaussSums[1]/GaussSums[0]]
    
valuationsPos1 = []
literalDifferencePos1 = []
for i in range(0, len(RatioComparisonsPos1)):
    diff = RatioComparisonsPos1[i] - GaussSumRatios[i]
    literalDifferencePos1.append(diff)
    val = (diff).valuation()
    valuationsPos1.append(val)

print "valuationsPos1"
print  valuationsPos1
    


