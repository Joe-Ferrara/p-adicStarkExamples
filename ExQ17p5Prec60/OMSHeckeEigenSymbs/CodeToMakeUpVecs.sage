## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

p = 5
for i in range(0, 32):
    phi = load('phi' + str(i))
    phiUp = phi.hecke(p)
    save(phiUp, 'phi' + str(i) + 'Up')