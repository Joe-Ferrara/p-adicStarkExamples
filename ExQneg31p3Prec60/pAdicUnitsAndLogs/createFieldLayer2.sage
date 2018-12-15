C.<zeta27> = CyclotomicField(9)
g = (zeta27 - 1).minpoly()

p = 3
prec = 60
ramInd = 18
sizeResField = 27
pAdics = Qp(p, prec)
Rp.<X> = PolynomialRing(pAdics)
g = Rp(g)
CpUnif.<D> = pAdics.extension(g)
Zeta9 = D + 1 ## THIS LINE IS A CHANGE TO THE OLD CODE

d = -31
logN = 2*prec*ramInd + 1
## THIS IS THE NUMBER OF TERMS IN THE SUM DEFINING LOG THAT WILL BE SUMMED
## I THINK 2*prec + 1 IS THE RIGHT NUMBER OF TERMS TO SUM BUT I'M NOT SURE

class Mp:
    def __init__(self, a, b):
        self._a = CpUnif(a); self._b = CpUnif(b)

    def __repr__(self):
        return '%s + (%s)*Alpha' %(self._a, self._b)

    def __add__(self, right):
        return Mp(self._a + right._a, self._b + right._b)

    def __sub__(self, right):
        return Mp(self._a - right._a, self._b - right._b)

    def __neg__(self):
        return Mp(-self._a, -self._b)

    def __mul__(self, right):
        aProd = self._a*right._a + d*self._b*right._b
        bProd = self._a*right._b + self._b*right._a
        return Mp(aProd, bProd)

    def __div__(self, right):
        norm = right._a*right._a - d*right._b*right._b
        numerator1 = self._a*right._a - d*self._b*right._b
        numerator2 = self._b*right._a - self._a*right._b
        return Mp(numerator1/norm, numerator2/norm)

    def __pow__(self, exponent):
        ans = Mp(1, 0)
        for i in range(0, exponent):
            ans = ans*self
        return ans

    def firstComponent(self):
        return self._a

    def secondComponent(self):
        return self._b

    def log_p(self):
        oneUnit = self**(sizeResField - 1)
        negOne = Mp(-1, 0)
        logOneUnit = Mp(0,0)
        variable = oneUnit + negOne
        for i in range(0, logN):
            logOneUnit = (negOne**i)*(variable**(i + 1))/Mp(i + 1, 0) + logOneUnit
        ans = logOneUnit/Mp(sizeResField - 1, 0)
        return ans

    def project(self, sign):
        ubar = Mp(self._a, CpUnif(-self._b))
        if sign == 1:
            return u*ubar
        if sign == -1:
            return u/ubar

