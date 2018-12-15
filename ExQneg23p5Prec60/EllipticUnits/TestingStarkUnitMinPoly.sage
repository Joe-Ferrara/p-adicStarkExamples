R.<x> = PolynomialRing(QQ)
F.<a> = NumberField(x^2 + 23)
S.<y> = PolynomialRing(F)
H.<b> = F.extension(y^3 - y + 1)
T.<z> = PolynomialRing(H)
H1.<c> = H.extension(z^5 - 10*z^3 + 5*z^2 + 10*z + 1)
U.<w> = PolynomialRing(H1)


f = w^15 - 832535*w^14 + 65231675*w^13 - 5650639400*w^12 + 15533478425*w^11 - 39376942640*w^10 - 212804236525*w^9 - 380541320125*w^8 - 2607229594750*w^7 - 2183192838625*w^6 + 3771011381950*w^5 - 1207366794625*w^4 + 99067277500*w^3 - 221569375*w^2 + 466875*w - 125
