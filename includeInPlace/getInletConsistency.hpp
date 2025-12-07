// Known distributions must be read above
// Reads prescribed floats ux, uy, uz from local cell
// Calculates float rho from consistency condition
// Works for face cells only. Outer normal points in the direction where no neighbours are

ux = uxArrayView[cell];
uy = uyArrayView[cell];
uz = uzArrayView[cell];

if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + f3 + f4 + f5 + f6 + f13 + f14 + f17 + f18;
    const float rho = fkProduct / (1.f - ux);
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + f1 + f2 + f3 + f4 + f7 + f8 + f9 + f10;
    const float rho = fkProduct / (1.f - uy);
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
{
	const float fkProduct = + f0 + f1 + f2 + f5 + f6 + f11 + f12 + f15 + f16;
    const float rho = fkProduct / (1.f - uz);
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + 2.f * f2 + f3 + f4 + f5 + f6 + 2.f * f8 + 2.f * f10 + 2.f * f11 + f13 + f14 + 2.f * f15 + f17 + f18 + 2.f * f19 + 2.f * f21 + 2.f * f24 + 2.f * f25;
    const float rho = fkProduct / (1.f - ux);
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + f1 + f2 + f3 + f4 + 2.f * f5 + f7 + f8 + f9 + f10 + 2.f * f11 + 2.f * f14 + 2.f * f16 + 2.f * f18 + 2.f * f20 + 2.f * f21 + 2.f * f23 + 2.f * f25;
    const float rho = fkProduct / (1.f - uy);
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
{
	const float fkProduct = + f0 + f1 + f2 + 2.f * f3 + f5 + f6 + 2.f * f7 + 2.f * f10 + f11 + f12 + 2.f * f13 + f15 + f16 + 2.f * f18 + 2.f * f19 + 2.f * f22 + 2.f * f23 + 2.f * f25;
    const float rho = fkProduct / (1.f - uz);
}
