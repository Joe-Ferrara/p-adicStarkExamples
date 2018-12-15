## the characteristic polynomial of T3 acting on the U5 eigenspace with eigenvalue 1 (resp -1) is x^2 + 2x (resp x^2 - 2x). Therefore is phi has U5 eigenvalue 1 (resp -1), then phi|T3 - 2phi (resp phi|T3 + 2phi) has T3 eigenvalue 0. 

upEigenSymbs = load("upMinusEigenSymbs")

heckeEigenSymbs = []
phi = upEigenSymbs[0]
eta = phi.hecke(3) + phi.scale(2)
heckeEigenSymbs.append(eta)
phi = upEigenSymbs[1]
eta = phi.hecke(3) + phi.scale(-2)
heckeEigenSymbs.append(eta)
save(heckeEigenSymbs, 'minusHeckeEigenSymbs')

checks = []
for q in prime_range(2, 6):
    temp = [q]
    for phi in heckeEigenSymbs:
        temp.append(phi.is_Tq_eigen(q))
    checks.append(temp)
save(checks, 'checkHeckeEigenSymbs')
