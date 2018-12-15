def make_St_units_from_gRoot(eta):
    u = (eta + (eta^2 - 4).sqrt())/2
    ubar = (eta - (eta^2 - 4).sqrt())/2
    return [u, ubar]

StarkUnits = []
for root in gRoots:
    eta = root[0]
    units = make_St_units_from_gRoot(eta)
    StarkUnits.append(units[0])
    StarkUnits.append(units[1])

## save stuff

save(StarkUnits, 'algebraicStarkUnits')