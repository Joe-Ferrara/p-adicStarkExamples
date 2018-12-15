precisions = load('newPrecisions')
N = min(precisions)
N

basis = []
for i in range(0, 32):
    phi = load('phi' + str(i))
    phi = phi.change_precision(N)
    basis.append(phi)

A = Matrix(Zmod(5), [phi.vector_of_total_measures() for phi in basis])
A.rank()

save(basis, 'minusBasisNewPrec')
