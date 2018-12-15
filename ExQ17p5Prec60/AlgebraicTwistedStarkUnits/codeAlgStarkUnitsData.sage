
## load

StarkUnits = load("orderedStarkUnits")

## make lists of the coefficients of the Stark Units

orderedAlgStarkUnitsData = []

for i in range(0, len(StarkUnits)):
    print i
    u = StarkUnits[i]
    temp = []
    coeffsCycUnif = u.list()
    for i in range(0, len(coeffsCycUnif)):
        coeffsQ = coeffsCycUnif[i].list()
        temp.append(coeffsQ)
    orderedAlgStarkUnitsData.append(temp)

## save stuff

save(orderedAlgStarkUnitsData, 'orderedAlgStarkUnitsData')