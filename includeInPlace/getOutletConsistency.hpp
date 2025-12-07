// Known distributions must be read above
// Reads prescribed float rho from local cell
// Calculates normal float u_ from consistency condition
// Calculates two tangential u_, u_ from moments
// Works for face cells only. Outer normal points in the direction where no neighbours are

const float rho = rhoArrayView[cell];

if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + f3 + f4 + f5 + f6 + f13 + f14 + f17 + f18;
    float ux = 1.f - fkProduct / rho;
	const float fkTangential1 = - f5 + f6 + f13 - f14 + f17 - f18;
	const float fkTangential2 = - f3 + f4 - f13 + f14 + f17 - f18;
    float uy = fkTangential1 / ( (2.f/3.f - ux * ux ) * rho );
    float uz = fkTangential2 / ( (2.f/3.f - ux * ux ) * rho );
}
elif (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + f1 + f2 + f3 + f4 + f7 + f8 + f9 + f10;
    float uy = 1.f - fkProduct / rho;
	const float fkTangential1 = + f1 - f2 + f7 - f8 + f9 - f10;
	const float fkTangential2 = - f3 + f4 - f7 + f8 + f9 - f10;
    float ux = fkTangential1 / ( (2.f/3.f - uy * uy ) * rho );
    float uz = fkTangential2 / ( (2.f/3.f - uy * uy ) * rho );
}
elif (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
{
	const float fkProduct = + f0 + f1 + f2 + f5 + f6 + f11 + f12 + f15 + f16;
    float uz = 1.f - fkProduct / rho;
	const float fkTangential1 = + f1 - f2 - f11 + f12 - f15 + f16;
	const float fkTangential2 = - f5 + f6 - f11 + f12 + f15 - f16;
    float ux = fkTangential1 / ( (2.f/3.f - uz * uz ) * rho );
    float uy = fkTangential2 / ( (2.f/3.f - uz * uz ) * rho );
}
elif (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + 2.f * f2 + f3 + f4 + f5 + f6 + 2.f * f8 + 2.f * f10 + 2.f * f11 + f13 + f14 + 2.f * f15 + f17 + f18 + 2.f * f19 + 2.f * f21 + 2.f * f24 + 2.f * f25;
    float ux = 1.f - fkProduct / rho;
	const float fkTangential1 = - f5 + f6 + f13 - f14 + f17 - f18;
	const float fkTangential2 = - f3 + f4 - f13 + f14 + f17 - f18;
    float uy = fkTangential1 / ( (2.f/3.f - ux * ux ) * rho );
    float uz = fkTangential2 / ( (2.f/3.f - ux * ux ) * rho );
}
elif (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
{
	const float fkProduct = + f0 + f1 + f2 + f3 + f4 + 2.f * f5 + f7 + f8 + f9 + f10 + 2.f * f11 + 2.f * f14 + 2.f * f16 + 2.f * f18 + 2.f * f20 + 2.f * f21 + 2.f * f23 + 2.f * f25;
    float uy = 1.f - fkProduct / rho;
	const float fkTangential1 = + f1 - f2 + f7 - f8 + f9 - f10;
	const float fkTangential2 = - f3 + f4 - f7 + f8 + f9 - f10;
    float ux = fkTangential1 / ( (2.f/3.f - uy * uy ) * rho );
    float uz = fkTangential2 / ( (2.f/3.f - uy * uy ) * rho );
}
elif (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
{
	const float fkProduct = + f0 + f1 + f2 + 2.f * f3 + f5 + f6 + 2.f * f7 + 2.f * f10 + f11 + f12 + 2.f * f13 + f15 + f16 + 2.f * f18 + 2.f * f19 + 2.f * f22 + 2.f * f23 + 2.f * f25;
    float uz = 1.f - fkProduct / rho;
	const float fkTangential1 = + f1 - f2 - f11 + f12 - f15 + f16;
	const float fkTangential2 = - f5 + f6 - f11 + f12 + f15 - f16;
    float ux = fkTangential1 / ( (2.f/3.f - uz * uz ) * rho );
    float uy = fkTangential2 / ( (2.f/3.f - uz * uz ) * rho );
}
