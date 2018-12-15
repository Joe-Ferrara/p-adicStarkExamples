



C.<zeta25> = CyclotomicField(25)
d = zeta25 - 1
prec = 60
ramInd = 20
p = 5

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




ratioComps1Pos1 = load("RatioComparisonsPos1NoGauss1")
ratioComps1Neg1 = load("RatioComparisonsNeg1NoGauss1")
approx1Pos1 = []
for ratioPair in ratioComps1Pos1:
    ratio = ratioPair[0]
    data = [ratio.valuation(), ratio.list()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Pos1.append(approx)
approx1Neg1 = []
for ratioPair in ratioComps1Neg1:
    ratio = ratioPair[0]
    data = [ratio.valuation(), ratio.list()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx1Neg1.append(approx)


ratioComps2Pos1 = load("RatioComparisonsPos1NoGauss2")
ratioComps2Neg1 = load("RatioComparisonsNeg1NoGauss2")
approx2Pos1 = []
for ratioPair in ratioComps2Pos1:
    ratio = ratioPair[0]
    data = [ratio.valuation(), ratio.list()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx2Pos1.append(approx)
approx2Neg1 = []
for ratioPair in ratioComps2Neg1:
    ratio = ratioPair[0]
    data = [ratio.valuation(), ratio.list()]
    alg = make_algebraic(data)
    approx = approximate(alg)
    approx2Neg1.append(approx)







