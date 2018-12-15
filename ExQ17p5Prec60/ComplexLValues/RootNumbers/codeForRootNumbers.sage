from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding

load("Functions.sage")

########################################################################################
##  Words about "notation" and strategy. The local root number at a place v such that ##
##  v^e divides the conuctor of the character with e > 0, is given by the sum:        ##
##      chi_v(beta)/IdealNorm(v^(e/2)) * (multiply by sum in next line)               ##
##             sum_(x in (O_(K_v)/v^e)^*) chi_v^(-1)(x)exp(2*pi*i LocalTrace(x/beta)) ##
##  where chi_v is the local at v part of the idelic character of chi                 ##
##        beta is an element of v^e - v^(e + 1)                                       ##
########################################################################################
##  Second note is that while the overarching strategy and outline is general, all    ##
##  of the code is hardwired to the specific example given in the parameters section  ##
########################################################################################
##  Third note: the conuctor of the character(s) we consider is                       ##
##          cond = (25)*((1/2)*alpha + 1/2) = (5)^2*((-1/2)*alpha + 3/2)              ##
##  where alpha is a square root of 17. These numbers are hardwired into the code at  ##
##  times.                                                                            ##
########################################################################################

def make_beta(K, L, prime, cond):
    ## beta in this case is a generator
    ## for the prime powers that divide
    ## the conductor
    if prime == K.ideal(5):
        return 25
    else:
        alpha = K.gen()
        return (-3/2)*alpha + 13/2
        ## NOTE: We use the generator
        ## (-3/2)*a + 13/2 = ((-1/2)*alpha + 3/2)*(alpha - 4)
        ## (alpha - 4 is a fundamental unit)
        ## because it has norm 4 making a later calculation
        ## easier.
        ## (-1/2)*alpha + 3/2 has norm -4

def sgn(element, K):
    ## K is a quadratic field
    ## this determines the sign
    ## (plus or minus 1)
    ## of the element of K using the
    ## negative suare root of the
    ## discriminant
    ## NOTE: We use the negative
    ## squareroot because the sign
    ## needs to be checked at the
    ## place that ramifies.
    d = K.discriminant()
    alpha = d.sqrt().n()
    sig = NumberFieldEmbedding(K, RLF, -alpha)
    if sig(element) > 0:
        return 1
    if sig(element) < 0:
        return -1

def chi_beta(K, L, psi, prime, beta):
    ## determines the value of the local character of chi at prime evaluated on beta
    ## the strategy is to use that chi as an idelic character is trivial on elements
    ## of K and in this example both betas are elements of K (so global elements)
    if prime == K.ideal(5):
        return 1
        ## Two notes:
        ## One: 5 is a uniformizer at K.ideal(5),
        ##      so chi_5(25) will not contribute
        ##      anything at any unramified place
        ## Two: 5 (and so 25) is 1 mod the other
        ##      prime power dividing the conductor
        ##      therefore chi_5(25) doesn't contribute
        ##      anything at the other prime dividing
        ##      the conductor
    else:
        ## in this case beta still generates the
        ## ideal prime, as above, but beta is
        ## not 1 mod 25, so we need to be clever
        ## to evaluate the local character at beta
        reslist = [beta, 1]
        Ilist = [K.ideal(25), K.ideal(beta)]
        gamma = K.solve_CRT(reslist, Ilist, True)
        epsilon = beta*gamma
        ideal = K.ideal(gamma)
        ans = psi(ideal.absolute_norm())*chiArtin(K, L, ideal, 1)*sgn(epsilon, K)
        return ans


def correct_invertible_residues(K, L, prime, cond):
    ## make a list of invertible residues of prime^e (where prime^e
    ## exactly divides cond) such that the invertible residues are
    ## 1 mod the other prime power that divides cond
    e = cond.valuation(prime)
    invRes = list((prime^e).invertible_residues())
    if prime == K.ideal(5):
        ideal = L.relative_discriminant() ## the other prime power in cond
        ans = []
        for x in invRes:
            reslist = [x, 1]
            Ilist = [prime^e, ideal]
            y = K.solve_CRT(reslist, Ilist, True)
            ans.append(y)
        return ans
    else:
        ideal = K.ideal(25)
        ans = []
        for x in invRes:
            reslist = [x, 1]
            Ilist = [prime^e, ideal]
            y = K.solve_CRT(reslist, Ilist, True)
            ans.append(y)
        return ans

def chi_inv_res(K, L, psi, invRes):
    ## NOTE: this returnes a list of chi_v^(-1) of the invertible
    ##       residues at v (note that v is not one of the inputs
    ##       but v is implicit in invRes)
    ans = []
    for x in invRes:
        ideal = K.ideal(x)
        temp = chiArtin(K, L, ideal, 1)*psi(ideal.absolute_norm())*sgn(x, K)
        ans.append(temp)
    return ans

def roots_of_unity_in_sum(K, invRes, beta):
    if beta == 25:
        C.<zeta25> = CyclotomicField(25)
        ans = []
        for x in invRes:
            trace = x.trace()
            ans.append(zeta25^trace)
        return ans
    else:
        n = beta.norm() ## this n is (or should be) 4
        C.<zetan> = CyclotomicField(n)
        ans = []
        for x in invRes:
            exponent = x*(beta.galois_conjugate()) + (x.galois_conjugate())*beta
            ans.append(zetan^exponent)
        return ans


def root_local_ram(K, L, psi, prime, cond):
    beta = make_beta(K, L, prime, cond)
    chibeta = chi_beta(K, L, psi, prime, beta)
    invRes = correct_invertible_residues(K, L, prime, cond)
    chiInvRes = chi_inv_res(K, L, psi, invRes)
    roots_of_unity = roots_of_unity_in_sum(K, invRes, beta)
    sum_part = 0
    C.<zeta100> = CyclotomicField(100)
    for i in range(0, len(invRes)):
        sum_part = sum_part + C(chiInvRes[i])*C(roots_of_unity[i])
    constant = chibeta/(prime.norm()) ## both the e's here are 2
    ans = C(constant)*sum_part
    return ans



## parameters

R.<x> = PolynomialRing(QQ)
K.<alpha> = NumberField(x^2 - 17)
S.<y> = PolynomialRing(K)
L.<xi> = K.extension(y^2 - (4 + alpha))
discLK = L.relative_discriminant()
p = 5
cond = discLK*K.ideal(p^2)
########################################################
## note that the conductor is divisible by two primes ##
## and both the primes are principal                  ##
## these two facts make calculations easier           ##
########################################################
different = K.different()
C1.<zeta5> = CyclotomicField(5)
C2.<zeta4> = CyclotomicField(4)
C3.<zeta100> = CyclotomicField(100)
G = DirichletGroup(25, C1, zeta5, 5)

## we will make one root number for each element of G (G is a group of Dirichlet characters)

rootNumbers = []

for psi in G:

## the root number will be a list of local root numbers we multiply together

    Wterms = []

## root number at infinity is -i since the character is of mixed signature

    Winf = zeta4^3
    Wterms.append(C3(Winf))

## root number at places not dividing the conductor is the character evaluated on the different

    Wdifferent = psi(different.absolute_norm())*chiArtin(K, L, different, 1)
    Wterms.append(C3(Wdifferent))

## we now calculate the root number at each of the primes dividing the conductor

    prime_factors = cond.prime_factors()
    for prime in prime_factors:
        Wprime = root_local_ram(K, L, psi, prime, cond)
        Wterms.append(Wprime)
    W = prod(Wterms)
    rootNumbers.append(W)

save(rootNumbers, 'rootNumbers')
