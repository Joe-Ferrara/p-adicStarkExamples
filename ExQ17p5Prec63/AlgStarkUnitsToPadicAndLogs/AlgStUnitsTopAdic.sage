#####################################################
## this code takes the ordered algebraic stark     ##
## units data and produces the p-adic logarithms   ##
## of the stark units.                             ##
## computationally taking the logs takes some      ##
## time.                                           ##
## before taking p-adic log, we project the stark  ##
## units to the subspace corresponding to the      ##
## p-stabilization.                                ##
#####################################################
#####################################################



## load

load("createField.sage")
algStarkUnitsData = load("orderedAlgStarkUnitsData")


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

pAdicStarkUnits = []

for ulist in algStarkUnitsData:
    embedding = pAdicEmbeddings(ulist, Zeta25, A)
    pAdicStarkUnits.append(embedding)

print "done part 1"

## make the projectiongs

print "project"
projections = []
for u in pAdicStarkUnits:
    print "projecting positive 1"
    uProj = u.project(1)
    projections.append(uProj)

print "logarithms"
pAdicLogs = []
for i in range(0, len(projections)):
    print "logging positive 1"
    u = projections[i]
    logu = u.log_p()
    logu = logu.firstComponent()
    save(logu, "Pos1logu" + str(i))
    pAdicLogs.append(logu)
save(pAdicLogs, "pAdicLogsPos1Proj")

print "project"
projections = []
for u in pAdicStarkUnits:
    print "projecting negative 1"
    uProj = u.project(-1)
    projections.append(uProj)

print "logarithms"
pAdicLogs = []
for i in range(0, len(projections)):
    print "logging negative 1"
    u = projections[i]
    logu = u.log_p()
    logu = logu.secondComponent()
    save(logu, "Neg1logu" + str(i))
    pAdicLogs.append(logu)
save(pAdicLogs, "pAdicLogsNeg1Proj")


