
print 1

StLogsPos1Pre = load("pAdicLogsPos1Proj0")
StLogsNeg1Pre = load("pAdicLogsNeg1Proj0")

print 2

StLogsPos1Data = []
StLogsNeg1Data = []
for i in range(0, len(StLogsPos1Pre)):
    tempPos1 = StLogsPos1Pre[i]
    tempNeg1 = StLogsNeg1Pre[i]
    dataPos1 = [tempPos1.valuation(), tempPos1.list()]
    dataNeg1 = [tempNeg1.valuation(), tempNeg1.list()]
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

## make tau(psi^-k)*StSum(psi^k) which if in Q_p(zeta5) then there is a mistake

print 8

GSStSumsTestPos1 = []
GSStSumsTestNeg1 = []
for k in range(0, 4):
    GSStSumPos1 = GaussSums[3-k]*StSumsPos1[k]
    GSStSumsTestPos1.append(GSStSumPos1)
    GSStSumNeg1 = GaussSums[3-k]*StSumsNeg1[k]
    GSStSumsTestNeg1.append(GSStSumNeg1)



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

print 9

algPos1 = []
algNeg1 = []
algTestPos1 = []
algTestNeg1 = []
for k in range(0, 4):
    temp = GSStSumsPos1[k]
    data = [temp.valuation(), temp.list()]
    alg = make_algebraic(data)
    algPos1.append(alg)
    temp = GSStSumsNeg1[k]
    data = [temp.valuation(), temp.list()]
    alg = make_algebraic(data)
    algNeg1.append(alg)
    temp = GSStSumsTestPos1[k]
    data = [temp.valuation(), temp.list()]
    alg = make_algebraic(data)
    algTestPos1.append(alg)
    temp = GSStSumsTestNeg1[k]
    data = [temp.valuation(), temp.list()]
    alg = make_algebraic(data)
    algTestNeg1.append(alg)

convenience = [algPos1, algNeg1, algTestPos1, algTestNeg1]

print 10

convenientValuations = []
for i in range(0, 4):
    valuationsList = []
    algNums = convenience[i]
    for k in range(0, 4):
        alg = algNums[k]
        data = alg.list()
        valuations = []
        for j in range(0, len(data)):
            valuations.append(data[j].valuation(5))
        valuationsList.append(valuations)
    convenientValuations.append(valuationsList)

## these are the results tau(psi^k)*StSum(psi^k) are in the correct field and tau(psi^(-k))*StSum(psi^k) are not

## sage: convenientValuations[0]

## [[2, 60, 60, 61, 60, 2, 60, 60, 62, 60, 3, 60, 60, 61, 60, 2, 60, 60, 61, 60],
##  [2, 60, 60, 60, 60, 2, 60, 60, 60, 60, 2, 60, 60, 60, 60, 3, 60, 60, 60, 60],
##  [2, 60, 61, 60, 60, 2, 60, 61, 60, 60, 2, 60, 61, 60, 60, 2, 60, 61, 60, 60],
##  [2, 60, 60, 60, 61, 2, 60, 60, 60, 62, 2, 60, 60, 60, 61, 2, 60, 60, 60, 61]]
## sage: convenientValuations[1]

## [[3, 61, 61, 62, 61, 2, 61, 61, 61, 61, 2, 61, 62, 63, 61, 2, 61, 61, 61, 61],
##  [2, 61, 61, 62, 62, 2, 61, 61, 63, 61, 3, 61, 61, 61, 64, 2, 61, 61, 63, 61],
##  [2, 65, 61, 61, 62, 2, 61, 61, 61, 61, 2, 61, 61, 61, 61, 2, 62, 61, 62, 61],
##  [2, 61, 61, 61, 61, 2, 61, 61, 61, 61, 2, 61, 61, 61, 62, 3, 61, 61, 61, 61]]
## sage: convenientValuations[2]

## [[60, 2, 60, 60, 61, 60, 3, 60, 60, 61, 60, 2, 60, 60, 61, 60, 2, 60, 60, 61],
##  [60, 60, 2, 60, 60, 60, 60, 2, 60, 60, 60, 60, 2, 60, 60, 60, 60, 3, 60, 60],
##  [61, 60, 60, 2, 60, 61, 60, 60, 2, 60, 61, 60, 60, 2, 60, 61, 60, 60, 2, 60],
##  [60, 60, 60, 61, 2, 60, 60, 60, 61, 2, 60, 60, 60, 61, 2, 60, 60, 60, 61, 2]]
## sage: convenientValuations[3]

## [[61, 2, 64, 61, 62, 61, 2, 61, 61, 62, 61, 2, 61, 61, 62, 61, 3, 61, 62, 61],
##  [61, 61, 2, 61, 62, 61, 61, 2, 61, 63, 61, 61, 3, 61, 61, 61, 62, 2, 61, 62],
##  [61, 61, 61, 2, 61, 61, 61, 61, 3, 61, 61, 61, 61, 2, 62, 61, 61, 62, 2, 61],
##  [61, 62, 62, 61, 2, 61, 62, 61, 61, 2, 61, 61, 61, 61, 2, 61, 61, 61, 61, 2]]






