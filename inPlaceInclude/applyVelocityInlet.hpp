// Applies velocity inlet moment based boundary condition on a face cell
// 1. Reads known distributions and prescribed ux, uy, uz
// 2. Calculates rho from consistency condition
// 3. Applies general moment based boundary condition to find unknown distributions

// Source paper: Pavel Eichler, Radek Fucik, and Pavel Strachota.
// Investigation of mesoscopic boundary conditions for lattice boltzmann method in laminar flow problems.
// Computers & Mathematics with Applications, 173:87–101, 2024.

// Reading prescribed velocity inlet ux, uy, uz
float ux = uxArrayView[cell];
float uy = uyArrayView[cell];
float uz = uzArrayView[cell];
// Declaring rho
float rho = 1.f;

if (flag == 1655) // outer normal [1, 0, 0]
{
	// Reading known distributions fk
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f26 = f26ArrayView[shiftedIndex[26]];
	// Applying consistency condition to find rho
	const float fkProduct = + f0 + 2.f * f1 + f3 + f4 + f5 + f6 + 2.f * f7 + 2.f * f9 + 2.f * f12 + f13 + f14 + 2.f * f16 + f17 + f18 + 2.f * f20 + 2.f * f22 + 2.f * f23 + 2.f * f26;
    rho = fkProduct / (1.f + ux);
	// At this point rho, ux, uy, uz are known
	// Applying general MBBC
	// Multiply K fk
	const float kf0 = + f0 + f1 + f3 + f4 + f5 + f6 + f7 + f9 + f12 + f13 + f14 + f16 + f17 + f18 + f20 + f22 + f23 + f26;
	const float kf1 = - f5 + f6 + f12 + f13 - f14 - f16 + f17 - f18 - f20 + f22 - f23 + f26;
	const float kf2 = - f3 + f4 - f7 + f9 - f13 + f14 + f17 - f18 + f20 - f22 - f23 + f26;
	const float kf3 = + f5 + f6 + f12 + f13 + f14 + f16 + f17 + f18 + f20 + f22 + f23 + f26;
	const float kf4 = + f3 + f4 + f7 + f9 + f13 + f14 + f17 + f18 + f20 + f22 + f23 + f26;
	const float kf5 = - f13 - f14 + f17 + f18 - f20 - f22 + f23 + f26;
	const float kf6 = - f13 + f14 + f17 - f18 + f20 - f22 - f23 + f26;
	const float kf7 = + f13 - f14 + f17 - f18 - f20 + f22 - f23 + f26;
	const float kf8 = + f13 + f14 + f17 + f18 + f20 + f22 + f23 + f26;
	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * uy;
	const float m2 = rho * uz;
	const float m3 = (1.f/3.f) * rho + rho * uy * uy;
	const float m4 = (1.f/3.f) * rho + rho * uz * uz;
	const float m5 = rho * uy * uz;
	const float m6 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m7 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m8 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	// Subtract m - Kfk
	const float s0 = m0 - kf0;
	const float s1 = m1 - kf1;
	const float s2 = m2 - kf2;
	const float s3 = m3 - kf3;
	const float s4 = m4 - kf4;
	const float s5 = m5 - kf5;
	const float s6 = m6 - kf6;
	const float s7 = m7 - kf7;
	const float s8 = m8 - kf8;
	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	f2 = + s0 - s3 - s4 + s8;
	f8 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	f10 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	f11 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	f15 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	f19 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f21 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f24 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f25 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
}
else if (flag == 1565) // outer normal [0, 1, 0]
{
	// Reading known distributions fk
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f26 = f26ArrayView[shiftedIndex[26]];
	// Applying consistency condition to find rho
	const float fkProduct = + f0 + f1 + f2 + f3 + f4 + 2.f * f6 + f7 + f8 + f9 + f10 + 2.f * f12 + 2.f * f13 + 2.f * f15 + 2.f * f17 + 2.f * f19 + 2.f * f22 + 2.f * f24 + 2.f * f26;
    rho = fkProduct / (1.f + uy);
	// At this point rho, ux, uy, uz are known
	// Applying general MBBC
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f3 + f4 + f6 + f7 + f8 + f9 + f10 + f12 + f13 + f15 + f17 + f19 + f22 + f24 + f26;
	const float kf1 = + f1 - f2 + f7 - f8 + f9 - f10 + f12 - f15 - f19 + f22 - f24 + f26;
	const float kf2 = - f3 + f4 - f7 + f8 + f9 - f10 - f13 + f17 - f19 - f22 + f24 + f26;
	const float kf3 = + f1 + f2 + f7 + f8 + f9 + f10 + f12 + f15 + f19 + f22 + f24 + f26;
	const float kf4 = + f3 + f4 + f7 + f8 + f9 + f10 + f13 + f17 + f19 + f22 + f24 + f26;
	const float kf5 = - f7 - f8 + f9 + f10 + f19 - f22 - f24 + f26;
	const float kf6 = - f7 + f8 + f9 - f10 - f19 - f22 + f24 + f26;
	const float kf7 = + f7 - f8 + f9 - f10 - f19 + f22 - f24 + f26;
	const float kf8 = + f7 + f8 + f9 + f10 + f19 + f22 + f24 + f26;
	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uz;
	const float m3 = (1.f/3.f) * rho + rho * ux * ux;
	const float m4 = (1.f/3.f) * rho + rho * uz * uz;
	const float m5 = rho * ux * uz;
	const float m6 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m7 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m8 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	// Subtract m - Kfk
	const float s0 = m0 - kf0;
	const float s1 = m1 - kf1;
	const float s2 = m2 - kf2;
	const float s3 = m3 - kf3;
	const float s4 = m4 - kf4;
	const float s5 = m5 - kf5;
	const float s6 = m6 - kf6;
	const float s7 = m7 - kf7;
	const float s8 = m8 - kf8;
	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	f5 = + s0 - s3 - s4 + s8;
	f11 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	f14 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	f16 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	f18 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	f20 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f21 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f23 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f25 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
}
else if (flag == 1556) // outer normal [0, 0, 1]
{
	// Reading known distributions fk
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f26 = f26ArrayView[shiftedIndex[26]];
	// Applying consistency condition to find rho
	const float fkProduct = + f0 + f1 + f2 + 2.f * f4 + f5 + f6 + 2.f * f8 + 2.f * f9 + f11 + f12 + 2.f * f14 + f15 + f16 + 2.f * f17 + 2.f * f20 + 2.f * f21 + 2.f * f24 + 2.f * f26;
    rho = fkProduct / (1.f + uz);
	// At this point rho, ux, uy, uz are known
	// Applying general MBBC
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f4 + f5 + f6 + f8 + f9 + f11 + f12 + f14 + f15 + f16 + f17 + f20 + f21 + f24 + f26;
	const float kf1 = + f1 - f2 - f8 + f9 - f11 + f12 - f15 + f16 + f20 - f21 - f24 + f26;
	const float kf2 = - f5 + f6 - f11 + f12 - f14 + f15 - f16 + f17 - f20 - f21 + f24 + f26;
	const float kf3 = + f1 + f2 + f8 + f9 + f11 + f12 + f15 + f16 + f20 + f21 + f24 + f26;
	const float kf4 = + f5 + f6 + f11 + f12 + f14 + f15 + f16 + f17 + f20 + f21 + f24 + f26;
	const float kf5 = + f11 + f12 - f15 - f16 - f20 + f21 - f24 + f26;
	const float kf6 = - f11 + f12 + f15 - f16 - f20 - f21 + f24 + f26;
	const float kf7 = - f11 + f12 - f15 + f16 + f20 - f21 - f24 + f26;
	const float kf8 = + f11 + f12 + f15 + f16 + f20 + f21 + f24 + f26;
	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = (1.f/3.f) * rho + rho * ux * ux;
	const float m4 = (1.f/3.f) * rho + rho * uy * uy;
	const float m5 = rho * ux * uy;
	const float m6 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m7 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m8 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;
	// Subtract m - Kfk
	const float s0 = m0 - kf0;
	const float s1 = m1 - kf1;
	const float s2 = m2 - kf2;
	const float s3 = m3 - kf3;
	const float s4 = m4 - kf4;
	const float s5 = m5 - kf5;
	const float s6 = m6 - kf6;
	const float s7 = m7 - kf7;
	const float s8 = m8 - kf8;
	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	f3 = + s0 - s3 - s4 + s8;
	f7 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	f10 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	f13 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	f18 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	f19 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f22 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f23 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f25 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
}
else if (flag == 1455) // outer normal [-1, 0, 0]
{
	// Reading known distributions fk
	f0 = f0ArrayView[shiftedIndex[0]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f17 = f17ArrayView[shiftedIndex[17]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f24 = f24ArrayView[shiftedIndex[24]];
	f25 = f25ArrayView[shiftedIndex[25]];
	// Applying consistency condition to find rho
	const float fkProduct = + f0 + 2.f * f2 + f3 + f4 + f5 + f6 + 2.f * f8 + 2.f * f10 + 2.f * f11 + f13 + f14 + 2.f * f15 + f17 + f18 + 2.f * f19 + 2.f * f21 + 2.f * f24 + 2.f * f25;
    rho = fkProduct / (1.f - ux);
	// At this point rho, ux, uy, uz are known
	// Applying general MBBC
	// Multiply K fk
	const float kf0 = + f0 + f2 + f3 + f4 + f5 + f6 + f8 + f10 + f11 + f13 + f14 + f15 + f17 + f18 + f19 + f21 + f24 + f25;
	const float kf1 = - f5 + f6 - f11 + f13 - f14 + f15 + f17 - f18 + f19 - f21 + f24 - f25;
	const float kf2 = - f3 + f4 + f8 - f10 - f13 + f14 + f17 - f18 - f19 + f21 + f24 - f25;
	const float kf3 = + f5 + f6 + f11 + f13 + f14 + f15 + f17 + f18 + f19 + f21 + f24 + f25;
	const float kf4 = + f3 + f4 + f8 + f10 + f13 + f14 + f17 + f18 + f19 + f21 + f24 + f25;
	const float kf5 = - f13 - f14 + f17 + f18 - f19 - f21 + f24 + f25;
	const float kf6 = - f13 + f14 + f17 - f18 - f19 + f21 + f24 - f25;
	const float kf7 = + f13 - f14 + f17 - f18 + f19 - f21 + f24 - f25;
	const float kf8 = + f13 + f14 + f17 + f18 + f19 + f21 + f24 + f25;
	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * uy;
	const float m2 = rho * uz;
	const float m3 = (1.f/3.f) * rho + rho * uy * uy;
	const float m4 = (1.f/3.f) * rho + rho * uz * uz;
	const float m5 = rho * uy * uz;
	const float m6 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m7 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m8 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	// Subtract m - Kfk
	const float s0 = m0 - kf0;
	const float s1 = m1 - kf1;
	const float s2 = m2 - kf2;
	const float s3 = m3 - kf3;
	const float s4 = m4 - kf4;
	const float s5 = m5 - kf5;
	const float s6 = m6 - kf6;
	const float s7 = m7 - kf7;
	const float s8 = m8 - kf8;
	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	f1 = + s0 - s3 - s4 + s8;
	f7 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	f9 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	f12 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	f16 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	f20 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f22 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f23 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f26 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
}
else if (flag == 1545) // outer normal [0, -1, 0]
{
	// Reading known distributions fk
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f4 = f4ArrayView[shiftedIndex[4]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f8 = f8ArrayView[shiftedIndex[8]];
	f9 = f9ArrayView[shiftedIndex[9]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f14 = f14ArrayView[shiftedIndex[14]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f20 = f20ArrayView[shiftedIndex[20]];
	f21 = f21ArrayView[shiftedIndex[21]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f25 = f25ArrayView[shiftedIndex[25]];
	// Applying consistency condition to find rho
	const float fkProduct = + f0 + f1 + f2 + f3 + f4 + 2.f * f5 + f7 + f8 + f9 + f10 + 2.f * f11 + 2.f * f14 + 2.f * f16 + 2.f * f18 + 2.f * f20 + 2.f * f21 + 2.f * f23 + 2.f * f25;
    rho = fkProduct / (1.f - uy);
	// At this point rho, ux, uy, uz are known
	// Applying general MBBC
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f3 + f4 + f5 + f7 + f8 + f9 + f10 + f11 + f14 + f16 + f18 + f20 + f21 + f23 + f25;
	const float kf1 = + f1 - f2 + f7 - f8 + f9 - f10 - f11 + f16 + f20 - f21 + f23 - f25;
	const float kf2 = - f3 + f4 - f7 + f8 + f9 - f10 + f14 - f18 + f20 + f21 - f23 - f25;
	const float kf3 = + f1 + f2 + f7 + f8 + f9 + f10 + f11 + f16 + f20 + f21 + f23 + f25;
	const float kf4 = + f3 + f4 + f7 + f8 + f9 + f10 + f14 + f18 + f20 + f21 + f23 + f25;
	const float kf5 = - f7 - f8 + f9 + f10 + f20 - f21 - f23 + f25;
	const float kf6 = - f7 + f8 + f9 - f10 + f20 + f21 - f23 - f25;
	const float kf7 = + f7 - f8 + f9 - f10 + f20 - f21 + f23 - f25;
	const float kf8 = + f7 + f8 + f9 + f10 + f20 + f21 + f23 + f25;
	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uz;
	const float m3 = (1.f/3.f) * rho + rho * ux * ux;
	const float m4 = (1.f/3.f) * rho + rho * uz * uz;
	const float m5 = rho * ux * uz;
	const float m6 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m7 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m8 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	// Subtract m - Kfk
	const float s0 = m0 - kf0;
	const float s1 = m1 - kf1;
	const float s2 = m2 - kf2;
	const float s3 = m3 - kf3;
	const float s4 = m4 - kf4;
	const float s5 = m5 - kf5;
	const float s6 = m6 - kf6;
	const float s7 = m7 - kf7;
	const float s8 = m8 - kf8;
	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	f6 = + s0 - s3 - s4 + s8;
	f12 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	f13 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	f15 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	f17 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	f19 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f22 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f24 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f26 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
}
else if (flag == 1554) // outer normal [0, 0, -1]
{
	// Reading known distributions fk
	f0 = f0ArrayView[shiftedIndex[0]];
	f1 = f1ArrayView[shiftedIndex[1]];
	f2 = f2ArrayView[shiftedIndex[2]];
	f3 = f3ArrayView[shiftedIndex[3]];
	f5 = f5ArrayView[shiftedIndex[5]];
	f6 = f6ArrayView[shiftedIndex[6]];
	f7 = f7ArrayView[shiftedIndex[7]];
	f10 = f10ArrayView[shiftedIndex[10]];
	f11 = f11ArrayView[shiftedIndex[11]];
	f12 = f12ArrayView[shiftedIndex[12]];
	f13 = f13ArrayView[shiftedIndex[13]];
	f15 = f15ArrayView[shiftedIndex[15]];
	f16 = f16ArrayView[shiftedIndex[16]];
	f18 = f18ArrayView[shiftedIndex[18]];
	f19 = f19ArrayView[shiftedIndex[19]];
	f22 = f22ArrayView[shiftedIndex[22]];
	f23 = f23ArrayView[shiftedIndex[23]];
	f25 = f25ArrayView[shiftedIndex[25]];
	// Applying consistency condition to find rho
	const float fkProduct = + f0 + f1 + f2 + 2.f * f3 + f5 + f6 + 2.f * f7 + 2.f * f10 + f11 + f12 + 2.f * f13 + f15 + f16 + 2.f * f18 + 2.f * f19 + 2.f * f22 + 2.f * f23 + 2.f * f25;
    rho = fkProduct / (1.f - uz);
	// At this point rho, ux, uy, uz are known
	// Applying general MBBC
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f3 + f5 + f6 + f7 + f10 + f11 + f12 + f13 + f15 + f16 + f18 + f19 + f22 + f23 + f25;
	const float kf1 = + f1 - f2 + f7 - f10 - f11 + f12 - f15 + f16 - f19 + f22 + f23 - f25;
	const float kf2 = - f5 + f6 - f11 + f12 + f13 + f15 - f16 - f18 + f19 + f22 - f23 - f25;
	const float kf3 = + f1 + f2 + f7 + f10 + f11 + f12 + f15 + f16 + f19 + f22 + f23 + f25;
	const float kf4 = + f5 + f6 + f11 + f12 + f13 + f15 + f16 + f18 + f19 + f22 + f23 + f25;
	const float kf5 = + f11 + f12 - f15 - f16 - f19 + f22 - f23 + f25;
	const float kf6 = - f11 + f12 + f15 - f16 + f19 + f22 - f23 - f25;
	const float kf7 = - f11 + f12 - f15 + f16 - f19 + f22 + f23 - f25;
	const float kf8 = + f11 + f12 + f15 + f16 + f19 + f22 + f23 + f25;
	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = (1.f/3.f) * rho + rho * ux * ux;
	const float m4 = (1.f/3.f) * rho + rho * uy * uy;
	const float m5 = rho * ux * uy;
	const float m6 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m7 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m8 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;
	// Subtract m - Kfk
	const float s0 = m0 - kf0;
	const float s1 = m1 - kf1;
	const float s2 = m2 - kf2;
	const float s3 = m3 - kf3;
	const float s4 = m4 - kf4;
	const float s5 = m5 - kf5;
	const float s6 = m6 - kf6;
	const float s7 = m7 - kf7;
	const float s8 = m8 - kf8;
	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	f4 = + s0 - s3 - s4 + s8;
	f8 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	f9 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	f14 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	f17 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	f20 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	f21 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f24 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	f26 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
}