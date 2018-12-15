## load

load("createField.sage")
algStarkUnitsData = load("orderedAlgStarkUnitsData")



###############################################
## NEED TO WRITE AN EXPLANATION OF THIS CODE ##
###############################################

def pAdicEmbeddings(ulist, zetaEmb, alphaEmb):
    ans = Mp(0,0)
    for i in range(0, len(ulist)):
        for j in range(0, len(ulist[i])):
            ans = Mp(ulist[i][j]*zetaEmb^j, 0)*alphaEmb^i + ans
    return ans


## make A and B which are (along with their negatives) the possible
## images of alpha and beta
## to determine a and b below we solve the equation
##    sqrt(4 + sqrt17) = const0 + const1*sqrt17
## for const0 and const1 in Qp

f = X^2 + 1
sqrtNeg1 = f.roots()[0][0]
## pick the square root of -1 that is 2 mod 5
if Integers(5)(sqrtNeg1) == Integers(5)(3):
    sqrtNeg1 = f.roots()[1][0]
f = X^2 - (2 - (1/2)*sqrtNeg1)
const0 = f.roots()[0][0]
## pick the root that is 1 mod 5
if Integers(5)(const0) == Integers(5)(4):
    const0 = f.roots()[1][0]
const1 = 1/(2*const0)
A = Mp(CpUnif(const0), CpUnif(const1))
B = Mp(CpUnif(const0), CpUnif(-const1))

## make the Stark units in all possible embeddings in Mp

##################################################################
## I only need the first 10 Stark units for theoretical reasons ##
## Therefore this code is different from that above             ##
##################################################################


## make all the embeddings

print "different embeddings"

alphaEmbeddings = [A, -A, B, -B]

pAdicStUnitsDiffEmb = []

for i in range(0, 5):
    print i
    zetaEmb = Zeta25^(2^(4*i))
    for j in range(0, 4):
        pAdicStarkUnitsL = []
        alphaEmb = alphaEmbeddings[j]
        for ulist in algStarkUnitsData:
            embedding = pAdicEmbeddings(ulist, zetaEmb, alphaEmb)
            pAdicStarkUnitsL.append(embedding)
        pAdicStUnitsDiffEmb.append(pAdicStarkUnitsL)

print "done part 1"

print len(pAdicStUnitsDiffEmb)


for i in range(0, len(pAdicStUnitsDiffEmb)):
    print i
    pAdicStarkUnits = pAdicStUnitsDiffEmb[i]

## make the projections

    print "project"
    projections = []
    for u in pAdicStarkUnits:
        uProj = u.project(1)
        projections.append(uProj)

    print "logarithms"
    pAdicLogs = []
    for u in projections:
        logu = u.log_p()
        logu = logu.firstComponent()
        pAdicLogs.append(logu)
    save(pAdicLogs, "pAdicLogsPos1Proj" + str(i))

    print "project"
    projections = []
    for u in pAdicStarkUnits:
        uProj = u.project(-1)
        projections.append(uProj)

    print "logarithms"
    pAdicLogs = []
    for u in projections:
        logu = u.log_p()
        logu = logu.secondComponent()
        pAdicLogs.append(logu)
    save(pAdicLogs, "pAdicLogsNeg1Proj" + str(i))


