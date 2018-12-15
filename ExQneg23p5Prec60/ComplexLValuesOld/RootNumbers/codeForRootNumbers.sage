
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

def make_beta(K, L, prime, cond):
    ## beta in this case is a generator
    ## for the prime powers that divide
    ## the conductor
    if prime == K.ideal(5):
        return 25

def chi_beta(K, L, psi, prime, beta):
    ## determines the value of the local character of chi at prime evaluated on beta
    ## the strategy is to use that chi as an idelic character is trivial on elements
    ## of K and in this example both betas are elements of K (so global elements)
    if prime == K.ideal(5):
        return 1
        ## 5 is a uniformizer at K.ideal(5),
        ## so chi_5(25) will not contribute
        ## anything at any unramified place

def chi_inv_res(K, L, psi, invRes):
    ## NOTE: this returnes a list of chi_v^(-1) of the invertible
    ##       residues at v (note that v is not one of the inputs
    ##       but v is implicit in invRes)
    ans = []
    C.<zeta75> = CyclotomicField(75)
    for x in invRes:
        ideal = K.ideal(x)
        temp = C(chiArtin(K, L, ideal, 1))*C(psi(ideal.absolute_norm()))
        ans.append(temp)
    return ans

def roots_of_unity_in_sum(K, invRes, beta):
    C.<zeta25> = CyclotomicField(25)
    ans = []
    for x in invRes:
        trace = x.trace()
        ans.append(zeta25^trace)
    return ans


def root_local_ram(K, L, psi, prime, cond):
    beta = make_beta(K, L, prime, cond)
    chibeta = chi_beta(K, L, psi, prime, beta)
    invRes = list(cond.invertible_residues())
    chiInvRes = chi_inv_res(K, L, psi, invRes)
    roots_of_unity = roots_of_unity_in_sum(K, invRes, beta)
    C.<zeta75> = CyclotomicField(75)
    sum_part = 0
    for i in range(0, len(invRes)):
        sum_part = sum_part + C(chiInvRes[i])*C(roots_of_unity[i])
    constant = chibeta/(prime.norm()) ## e is 2 here
    ans = constant*sum_part
    return ans



## parameters

R.<x> = PolynomialRing(QQ)
K.<alpha> = NumberField(x^2 + 23)
S.<y> = PolynomialRing(K)
L.<xi> = K.extension(y^3 - y^2 + 1)
p = 5
cond = K.ideal(p^2)
########################################################
## note that the conductor is divisible by one prime  ##
## which is principal                                 ##
## this makes calculations easier                     ##
########################################################
different = K.different()
C1.<zeta25> = CyclotomicField(25)
C2.<zeta75> = CyclotomicField(75)
zeta5 = zeta25^5
G = DirichletGroup(25, C, zeta5, 5)

## we will make one root number for element of G

rootNumbers = []

for psi in G:

## the root number will be a list of local root numbers we multiply together

    Wterms = []

## root number at infinity is 1 since the base field is imaginary quadratic

    Winf = 1
    Wterms.append(Winf)

## root number at places not dividing the conductor is the character evaluated on the different

    Wdifferent = C2(psi(different.absolute_norm()))*C2(chiArtin(K, L, different, 1))
    Wterms.append(C2(Wdifferent))

## we now calculate the root number at each of the primes dividing the conductor

    prime_factors = cond.prime_factors()
    for prime in prime_factors:
        Wprime = root_local_ram(K, L, psi, prime, cond)
        Wterms.append(Wprime)
    W = prod(Wterms)
    rootNumbers.append(W)

save(rootNumbers, 'rootNumbers')
