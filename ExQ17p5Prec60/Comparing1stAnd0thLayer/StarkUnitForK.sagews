︠55aa33a5-f118-44c9-8300-f5584a49b82a︠
R.<x> = PolynomialRing(QQ)
K.<b> = NumberField(x^4 - x^3 - 2*x^2 - x + 1, embedding=2.08101899662454)
RR.coerce_map_from(K)
︡3a53d9a7-b747-4383-a405-88b72b503ced︡{"stdout":"Composite map:\n  From: Number Field in b with defining polynomial x^4 - x^3 - 2*x^2 - x + 1\n  To:   Real Field with 53 bits of precision\n  Defn:   Generic morphism:\n          From: Number Field in b with defining polynomial x^4 - x^3 - 2*x^2 - x + 1\n          To:   Real Lazy Field\n          Defn: b -> 2.081018996624536?\n        then\n          Conversion via _mpfr_ method map:\n          From: Real Lazy Field\n          To:   Real Field with 53 bits of precision\n"}︡
︠1cbf3241-ba5f-4dd0-a921-bee39fc799bc︠
bnum = RR(b)
bnum
︡cf70785c-c9cf-4361-8bd4-8e3e2d875196︡{"stdout":"2.08101899662454\n"}︡
︠3c1da107-d0f0-4710-9801-9641875c2134︠
((1 + 17.sqrt())/4 +(((1 + 17.sqrt())/2).sqrt())/2).numerical_approx()
((1 + 17.sqrt())/2).numerical_approx()
︡8a5b0d71-8f58-4be6-a7ce-0cb1b0e5e8a1︡{"stdout":"2.08101899662454\n"}︡{"stdout":"2.56155281280883\n"}︡
︠68d8e613-1f6e-4e33-8133-12151fdacf5b︠
F.<a> = NumberField(x**2 - x - 4, embedding = 2.56155281280883)
RR.coerce_map_from(F)
︡630c218e-26cf-458a-95e4-9d1f067fe900︡{"stdout":"Composite map:\n  From: Number Field in a with defining polynomial x^2 - x - 4\n  To:   Real Field with 53 bits of precision\n  Defn:   Generic morphism:\n          From: Number Field in a with defining polynomial x^2 - x - 4\n          To:   Real Lazy Field\n          Defn: a -> 2.561552812808830?\n        then\n          Conversion via _mpfr_ method map:\n          From: Real Lazy Field\n          To:   Real Field with 53 bits of precision\n"}︡
︠2bdeb084-f5a0-4c1c-b32e-31d9830c97a3︠
RR(a)
︡4973a753-b857-4b3c-9e53-24fa7533580b︡{"stdout":"2.56155281280883\n"}︡
︠018c75b7-3b59-47fb-9f28-e83c4ad53f61︠
OK = K.OK()
OK
OK.basis()
OF = F.OK()
OF
OF.basis()
︡f8715654-d368-4c6d-9b84-b9920c70b36b︡{"stdout":"Maximal Order in Number Field in b with defining polynomial x^4 - x^3 - 2*x^2 - x + 1\n"}︡{"stdout":"[1, b, b^2, b^3]\n"}︡{"stdout":"Maximal Order in Number Field in a with defining polynomial x^2 - x - 4\n"}︡{"stdout":"[1, a]\n"}︡
︠d5efe220-2f40-4d8f-af26-27d16972fa42︠
OK(a)
︡3adce140-474c-4386-8cc7-95f5539760b2︡{"stdout":"-b^3 + b^2 + 3*b + 1"}︡{"stdout":"\n"}︡
︠348c043a-d944-4505-945c-5856578fa7ce︠
L1 = 0.5583995415
sqrt68num = 68.sqrt().numerical_approx()
pinum = pi.numerical_approx()
︡3a577228-6b05-4234-b7f0-50fe038d4df7︡
︠49c7e760-0308-483d-af73-ab7840c85bcc︠
Lderiv0 = (sqrt68num/pinum)*L1
Lderiv0
︡a0f6ea19-380d-4038-80d4-d37a11938fe1︡{"stdout":"1.46571535190609\n"}︡
︠f63e661a-698d-48c3-8799-d730930bbc11︠
UK = UnitGroup(K)
UK
UK.fundamental_units()
︡b6ded0dc-3d20-41e1-a13b-b1eab304fa80︡{"stdout":"Unit group with structure C2 x Z x Z of Number Field in b with defining polynomial x^4 - x^3 - 2*x^2 - x + 1\n"}︡{"stdout":"[b, b^3 - 2*b^2]\n"}︡
︠8c947b8f-5359-4d30-b307-443cdb656a0b︠
unum = RR(b)
vnum = RR(UK.fundamental_units()[1])
unum
vnum
︡5edd08a3-ffdb-47bb-b2bd-0be0023a2368︡{"stdout":"2.08101899662454\n"}︡{"stdout":"0.350864112752586\n"}︡
︠caa87738-49ea-44fc-a935-575b101f60c0︠
unum**4
(2*Lderiv0).exp()
︡a1425bed-e2c1-4e0f-b074-de0ecb7d8bc8︡{"stdout":"18.7544433666259\n"}︡{"stdout":"18.7544433650804\n"}︡
︠f7c21a95-9fd0-4787-86c9-a712a6fdda7a︠
Q5 = Qp(5, 12)
R5.<xx> = PolynomialRing(Q5)
F5.<aa> = Q5.extension(xx**2 - xx - 4)
F5
︡59af03fb-c376-4f81-bfa0-15d83117e903︡{"stdout":"Unramified Extension of 5-adic Field with capped relative precision 12 in aa defined by (1 + O(5^12))*xx^2 + (4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12))*xx + (1 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12))\n"}︡
︠f79150d3-7963-44fa-bcf7-1dac7ec177e7︠
aa
︡50078032-5e88-4349-a9b2-baed90a29994︡{"stdout":"aa + O(5^12)\n"}︡
︠c3fd1240-af4f-40ba-ab39-f6f280222b1c︠
I = K.ideal(5)
I.factor()
︡531e0b8f-ff8e-4c00-bba3-028d09c69e3c︡{"stdout":"(Fractional ideal (b^3 - 2*b^2 - b - 1)) * (Fractional ideal (-b^2 + 2*b + 2))\n"}︡
︠51ae6ba4-4142-4581-a208-0d176ef93e77︠
S0.<z> = GF(5^2, modulus = x^2 - x - 4)
︡1179b1e2-2eae-47b1-bbec-6a7895a6f585︡
︠d8cf5422-53ea-4eab-9528-6fc2ad634024︠
T.<X> = PolynomialRing(S0)
︡c8ff0c31-7e60-497b-9e22-c8c181cabd7c︡
︠d27a4d12-e948-4dab-9002-8f7cb967d484︠
f = X^2 - z*X + 1
f.factor()
︡8797b85d-d2be-4293-87b3-9d5f7f721cf7︡{"stdout":"(X + z + 1) * (X + 3*z + 4)\n"}︡
︠5085055d-cfb2-4290-9389-bee753397bdf︠
p = 5
S1.<z> = PolynomialRing(Zmod(p^2))
I = S1.ideal(z^2 - z - 4)
T1.<X> = PolynomialRing(S1.quotient_ring(I))
︡d4ecc92a-5a7e-46a2-9196-2361940796cf︡
︠34857f9f-93fc-41c9-989e-2ef37aa1f98c︠
f = X^2 - z*X + 1
︡0f1c3558-cf89-428b-8162-18ff086fa1c1︡
︠7775f1e0-cfb7-4dd8-bf34-ce6f5b4217c1︠
for i in range(0, p):
    for j in range(0, p):
        print [S1.quotient_ring(I)(1 + 2*z + i*p + j*z*p), f(S1.quotient_ring(I)(1 + 2*z + i*p + j*z*p))]
︡60ab7c8c-42d4-4332-8b4a-2bb0bbf59f33︡{"stdout":"[2*zbar + 1, 5*zbar + 10]\n[7*zbar + 1, 5*zbar + 20]\n[12*zbar + 1, 5*zbar + 5]\n[17*zbar + 1, 5*zbar + 15]\n[22*zbar + 1, 5*zbar]\n[2*zbar + 6, 20*zbar + 20]\n[7*zbar + 6, 20*zbar + 5]\n[12*zbar + 6, 20*zbar + 15]\n[17*zbar + 6, 20*zbar]\n[22*zbar + 6, 20*zbar + 10]\n[2*zbar + 11, 10*zbar + 5]\n[7*zbar + 11, 10*zbar + 15]\n[12*zbar + 11, 10*zbar]\n[17*zbar + 11, 10*zbar + 10]\n[22*zbar + 11, 10*zbar + 20]\n[2*zbar + 16, 15]\n[7*zbar + 16, 0]\n[12*zbar + 16, 10]\n[17*zbar + 16, 20]\n[22*zbar + 16, 5]\n[2*zbar + 21, 15*zbar]\n[7*zbar + 21, 15*zbar + 10]\n[12*zbar + 21, 15*zbar + 20]\n[17*zbar + 21, 15*zbar + 5]\n[22*zbar + 21, 15*zbar + 15]\n"}︡
︠79642be7-99a4-49ea-933a-a3928de21f50︠
S2.<z> = PolynomialRing(Zmod(p^3))
I = S2.ideal(z^2 - z - 4)
T2.<X> = PolynomialRing(S2.quotient_ring(I))
f = X^2 - z*X + 1
︡0d6b6fa0-7e61-4b09-a75c-9ef05823ec3f︡
︠56de884c-b833-4525-b7df-597baac2d7b8︠
for i in range(0, p):
    for j in range(0, p):
        print [S2.quotient_ring(I)(16 + 7*z + i*p^2 + j*z*p^2), f(S2.quotient_ring(I)(16 + 7*z + i*p^2 + j*z*p^2))]
︡5dca6490-b0e3-400e-9458-196bc1e23183︡{"stdout":"[7*zbar + 16, 50]\n[32*zbar + 16, 100]\n[57*zbar + 16, 25]\n[82*zbar + 16, 75]\n[107*zbar + 16, 0]\n[7*zbar + 41, 75*zbar + 100]\n[32*zbar + 41, 75*zbar + 25]\n[57*zbar + 41, 75*zbar + 75]\n[82*zbar + 41, 75*zbar]\n[107*zbar + 41, 75*zbar + 50]\n[7*zbar + 66, 25*zbar + 25]\n[32*zbar + 66, 25*zbar + 75]\n[57*zbar + 66, 25*zbar]\n[82*zbar + 66, 25*zbar + 50]\n[107*zbar + 66, 25*zbar + 100]\n[7*zbar + 91, 100*zbar + 75]\n[32*zbar + 91, 100*zbar]\n[57*zbar + 91, 100*zbar + 50]\n[82*zbar + 91, 100*zbar + 100]\n[107*zbar + 91, 100*zbar + 25]\n[7*zbar + 116, 50*zbar]\n[32*zbar + 116, 50*zbar + 50]\n[57*zbar + 116, 50*zbar + 100]\n[82*zbar + 116, 50*zbar + 25]\n[107*zbar + 116, 50*zbar + 75]\n"}︡
︠b67f0e59-5b59-4aba-97a4-ba08579f8478︠

︡39f0136c-3965-4b1b-b249-cc8a53709240︡
︠4cffde2f-e1b9-4e9d-b81e-82c279786470︠
RRR.<Z> = PolynomialRing(ZZ)
r = 1 + 2*Z
for k in range(2, 13):
    S.<z> = PolynomialRing(Zmod(p^k))
    I = S.ideal(z^2 - z - 4)
    T.<X> = PolynomialRing(S.quotient_ring(I))
    f = X^2 - z*X + 1
    for i in range(0, p):
        for j in range(0, p):
            t = r + i*p^(k-1) + j*Z*p^(k-1)
            s = S.quotient_ring(I)(r) + S.quotient_ring(I)(i*p^(k-1) + j*z*p^(k-1))
            if f(s) == 0:
                print s
                rr = t
    r = rr
︡7f897510-c055-4f6c-ab92-28d5aca3060d︡{"stdout":"7*zbar + 16\n107*zbar + 16"}︡{"stdout":"\n232*zbar + 516"}︡{"stdout":"\n2732*zbar + 3016"}︡{"stdout":"\n8982*zbar + 3016\n24607*zbar + 3016"}︡{"stdout":"\n337107*zbar + 3016"}︡{"stdout":"\n1118357*zbar + 1174891"}︡{"stdout":"\n8930857*zbar + 7034266"}︡{"stdout":"\n8930857*zbar + 26565516\n106587107*zbar + 124221766"}︡{"stdout":"\n"}︡
︠a5dd6c39-3699-4d20-b170-5bc69554907d︠










