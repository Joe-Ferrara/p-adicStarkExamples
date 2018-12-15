## FIRST MUST LOAD master.sage FROM DIFFERENT DIRECTORY

basis = load('plusBasisNewPrec')
upVecs = load('plusUpVecsNewPrec')

def makeUpMat(basis, upVecs):
    p = basis[0].p()
    prec = basis[0].num_moments()
    basisMat = Matrix(Zmod(p**prec), [phi.vector_of_total_measures() for phi in basis])
    matRows = []
    for n in range(0, len(basis)):
        b = vector(Zmod(p**(prec)), upVecs[n].vector_of_total_measures())
        x = basisMat.solve_left(b)
        matRows.append(x)
    ans = Matrix(Zmod(p**prec), [v for v in matRows])
    return ans

upMat = makeUpMat(basis, upVecs)
save(upMat, 'plusUpMat')