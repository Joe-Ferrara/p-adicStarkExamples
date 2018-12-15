## load

##########################################################
## NOTE WHAT THE PRECISION IS IN createFieldLayer1.sage ##
##########################################################

load("createFieldLayer1.sage")
algStarkUnitsData = load("orderedStarkUnitsData")

## the following function is hardwired for this example
## it works when ulist is generated as follows
## ulist = []
## for i in range(0, 5):
##     for j in range(0, 3):
##         for k in range(0, 2):
##             ulist.append(u.list()[i].list()[j].list()[k])
##


def pAdicEmbeddings(ulist, B, C):
    ans = Mp(0,0)
    for i in range(0, 5):
        for j in range(0, 3):
            ans = Mp(ulist[6*i + 2*j + 0], 0)*C^i*B^j + Mp(0, ulist[6*i + 2*j + 1])*C^i*B^j + ans
    return ans

## pick the root of y^3 - y + 1 that is in Q3
## call that root B

f = X^3 - X + 1
B = f.roots()[0][0]
B = Mp(CpUnif(B), 0)

## make C = Zeta25 + Zeta25^7 + Zeta25^(-1) + Zeta25^(-7)

C = Mp(Zeta25 + Zeta25^7 + Zeta25^(-1) + Zeta25^(-7), 0)

## make the embeddings of the Stark units

pAdicStUnits = []

for ulist in algStarkUnitsData:
    embedding = pAdicEmbeddings(ulist, B, C)
    pAdicStUnits.append(embedding)

## make the projections and logarithms
## Pos1 means projecting to space where complex conj acts as identity
## Neg1 means projecting to space where complex conj acts as inversion

projectionsPos1 = []
projectionsNeg1 = []
pAdicLogsPos1BeforeSimp = []
pAdicLogsNeg1BeforeSimp = []
pAdicLogsPos1 = []
pAdicLogsNeg1 = []

print "projecting"
for u in pAdicStUnits:
    print "poject"
    uProjPos1 = u.project(1)
    projectionsPos1.append(uProjPos1)
    uProjNeg1 = u.project(-1)
    projectionsNeg1.append(uProjNeg1)

print "logging"
for i in range(0, len(projectionsPos1)):
    print i
    u = projectionsPos1[i]
    logu = u.log_p()
    pAdicLogsPos1BeforeSimp.append(logu)
    logu = logu.firstComponent()
    pAdicLogsPos1.append(logu)
##    save(logu, "pAdicLogsPos1Proj" + str(i))
save(pAdicLogsPos1, "pAdicLogsPos1")

print "logging"
for i in range(0, len(projectionsNeg1)):
    print i
    u = projectionsNeg1[i]
    logu = u.log_p()
    pAdicLogsNeg1BeforeSimp.append(logu)
    logu = logu.secondComponent()
    pAdicLogsNeg1.append(logu)
##    save(logu, "pAdicLogsNeg1Proj" + str(i))
save(pAdicLogsNeg1, "pAdicLogsPos1")


pAdicLogsPos1Data = []
pAdicLogsNeg1Data = []
for i in range(0, len(pAdicLogsPos1)):
    tempPos1 = pAdicLogsPos1[i]
    tempNeg1 = pAdicLogsNeg1[i]
    dataPos1 = [tempPos1.valuation(), tempPos1.expansion()]
    dataNeg1 = [tempNeg1.valuation(), tempNeg1.expansion()]
    pAdicLogsPos1Data.append(dataPos1)
    pAdicLogsNeg1Data.append(dataNeg1)
save(pAdicLogsPos1Data, "pAdicLogsPos1Data")
save(pAdicLogsNeg1Data, "pAdicLogsNeg1Data")









 
