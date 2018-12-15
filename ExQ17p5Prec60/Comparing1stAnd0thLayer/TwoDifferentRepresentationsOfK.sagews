︠1b770437-601b-40ad-9ce5-6f0e3e71b1db︠

︡8b124e95-4d6b-4494-bb6e-30d3d579d556︡
︠b1e3dfc5-4ed6-4b2b-afbb-63e54eb64cf8︠
︡96c4a559-78ec-4046-8675-6c3acbef4262︡
︠b92797f0-8e99-403e-88ee-38a475afb4c4︠
pin = pi.numerical_approx(); pin
︡9539127a-82d2-4e1d-871d-2915aa2538c9︡{"stdout":"3.14159265358979\n"}︡
︠bce79a24-2671-44c8-adc8-ce77674caadf︠
nm = 68.sqrt().numerical_approx(); nm
︡716fc0b2-9eab-4d3c-bff0-62390370b69c︡{"stdout":"8.24621125123532\n"}︡
︠66a9c7fa-00a3-4f4d-a67d-e1402afb562a︠
L1 = 0.5583995415
︡87a7c6ef-7e9e-4e57-9e16-1558e4de3300︡
︠a1574603-afe3-4a6d-8a4a-19da8a3de7ac︠
L0 = (nm/pin)*L1
L0

︡26717861-7e24-4a54-954b-9110f46528ab︡{"stdout":"1.46571535190609\n"}︡
︠39426ee4-2e8c-4b03-8f66-bede04bb0b79︠
RR.<x> = PolynomialRing(QQ)
F.<a> = NumberField(x^2-17)
R.<y> = PolynomialRing(F)
K.<b> = F.extension(y^2-(4+a))
I = K.relative_discriminant()
I
NI = I.absolute_norm()
NI
D = F.discriminant()
D
N = D*NI
N
︡fa003202-bdfc-48de-920f-e39cb127a20d︡{"stdout":"Fractional ideal (1/2*a + 1/2)\n"}︡{"stdout":"4\n"}︡{"stdout":"17\n"}︡{"stdout":"68\n"}︡
︠7b145078-6da6-47c9-afdf-75ada072e016︠
FF.<A> = NumberField(x**2 - x - 4)
T.<Y> = PolynomialRing(FF)
KK.<B> = FF.extension(Y**2 - A*Y + 1)
S.<X> = PolynomialRing(K)
f = Y**4 - Y**3 - 2*Y**2 - Y + 1
f.factor()
g = X**4 - X**3 - 2*X**2 - X + 1
g.factor()
I = KK.relative_discriminant()
I
NI = I.absolute_norm()
NI
D = FF.discriminant()
D
N = D*NI
N
︡36a23fb4-a667-43e7-8afc-df45543f1127︡{"stdout":"(Y^2 - A*Y + 1) * (Y^2 + (A - 1)*Y + 1)\n"}︡{"stdout":"(X + (-1/4*a + 3/4)*b - 1/4*a - 1/4) * (X + (1/4*a - 3/4)*b - 1/4*a - 1/4) * (X^2 + (1/2*a - 1/2)*X + 1)\n"}︡{"stdout":"Fractional ideal (A)\n"}︡{"stdout":"4\n"}︡{"stdout":"17\n"}︡{"stdout":"68\n"}︡
︠08f11f83-0d42-41f7-a7c4-0bdde5294486︠
UK = UnitGroup(K); UK
UK.fundamental_units()

︡a647d776-b4ba-4ffe-b4c4-989b6bc82949︡{"stdout":"Unit group with structure C2 x Z x Z of Number Field in b with defining polynomial y^2 - a - 4 over its base field\n"}︡{"stdout":"[(1/4*a - 3/4)*b + 1/4*a + 1/4, b]\n"}︡
︠d294ca48-6db0-4bbc-bac0-954f731f66cf︠
aa = 17.sqrt().numerical_approx()
bb = (4 + aa).sqrt().numerical_approx()
︡85fc5def-eda4-4772-8c11-454a21db0460︡
︠f0e829a8-79cd-44fe-a833-3b6505ed9d53︠
vv = ((1/4*aa - 3/4)*bb + 1/4*aa + 1/4).numerical_approx()
uu = bb.numerical_approx()
uu
vv

︡3b497dd3-2b48-4db9-9440-0ef035ea7cf7︡{"stdout":"2.85010624812789\n"}︡{"stdout":"2.08101899662454\n"}︡
︠81dd68e6-ac6a-473a-9be9-46ccd182cfc2︠
aaa = uu.log()
aaa
bbb = vv.log()
bbb
ccc = 2*L0
ccc
︡29d01434-18ca-4321-8f3f-e3f2f40ff141︡{"stdout":"1.04735627363055\n"}︡{"stdout":"0.732857675973645\n"}︡{"stdout":"2.93143070381218\n"}︡
︠cc3568dd-9370-43b6-a203-cca049d56702︠
for x in range(0, 5):
    print x, ((ccc-aaa*x)/bbb).numerical_approx()
︡7bd71e17-0f00-4a98-b06a-c9940188e9bc︡{"stdout":"0 3.99999999988756\n1 2.57085992539891\n2 1.14171985091027\n3 -0.287420223578377\n4 -1.71656029806702\n"}︡
︠52ed2602-e472-4d98-89c9-afd5fa66cc3c︠
(2*L0).exp()
vv**4
︡1413e6a6-1430-4403-8871-b85fe27a73b2︡{"stdout":"18.7544433650804\n"}︡{"stdout":"18.7544433666259\n"}︡
︠b919d1d5-0d5e-4696-8158-97f2d002a8ea︠
v = (1/4*a - 3/4)*b + 1/4*a + 1/4
u = b
︡e4b8e338-5060-4388-ba25-bd4de451c90d︡
︠3af0fafe-08ba-42b1-a3e6-828486711a33︠
K.integral_basis()
︡f2089063-6e8a-4cff-89ba-d64e29fc85fe︡{"stdout":"[(1/4*a + 5/4)*b + 1/4*a + 5/4, (1/2*a + 5/2)*b, a + 4, (a + 4)*b]"}︡{"stdout":"\n"}︡
︠1a2a915a-5ee8-431f-b607-a1c3ec4acd2c︠










