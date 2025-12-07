// Known distribution floats f must be read above
// Floats rho, ux, uy, uz must be defined above
// Applies general moment based boundary condition on a boundary cell
// Calculates unknown distributions
// Outer normal identifies cell type - face cell / edge cell / corner cell and its orientation
// Example: cell has outer normal [1, -1, 0]
// -> Edge cell, has neighbours in both Z directions but no neighbours in +X and -Y direction

// Source paper: Pavel Eichler, Radek Fucik, and Pavel Strachota.
// Investigation of mesoscopic boundary conditions for lattice boltzmann method in laminar flow problems.
// Computers & Mathematics with Applications, 173:87–101, 2024.

if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0)
{
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
	float f2 = + s0 - s3 - s4 + s8;
	float f8 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	float f10 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	float f11 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	float f15 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	float f19 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f21 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f24 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f25 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0)
{
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
	float f5 = + s0 - s3 - s4 + s8;
	float f11 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	float f14 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	float f16 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	float f18 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	float f20 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f21 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f23 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f25 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1)
{
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
	float f3 = + s0 - s3 - s4 + s8;
	float f7 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	float f10 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	float f13 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	float f18 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	float f19 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f22 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f23 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f25 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0)
{
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
	float f1 = + s0 - s3 - s4 + s8;
	float f7 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	float f9 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	float f12 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	float f16 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	float f20 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f22 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f23 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f26 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0)
{
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
	float f6 = + s0 - s3 - s4 + s8;
	float f12 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	float f13 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	float f15 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	float f17 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	float f19 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f22 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f24 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f26 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
}
else if (outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1)
{
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
	float f4 = + s0 - s3 - s4 + s8;
	float f8 = - 0.5f * s1 + 0.5f * s3 + 0.5f * s7 - 0.5f * s8;
	float f9 = + 0.5f * s1 + 0.5f * s3 - 0.5f * s7 - 0.5f * s8;
	float f14 = - 0.5f * s2 + 0.5f * s4 + 0.5f * s6 - 0.5f * s8;
	float f17 = + 0.5f * s2 + 0.5f * s4 - 0.5f * s6 - 0.5f * s8;
	float f20 = - 0.25f * s5 - 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
	float f21 = + 0.25f * s5 - 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f24 = - 0.25f * s5 + 0.25f * s6 - 0.25f * s7 + 0.25f * s8;
	float f26 = + 0.25f * s5 + 0.25f * s6 + 0.25f * s7 + 0.25f * s8;
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f4 + f6 + f8 + f9 + f12 + f15 + f17 + f24 + f26;
	const float kf1 = + f1 - f2 - f8 + f9 + f12 - f15 - f24 + f26;
	const float kf2 = + f6 + f12 + f15 + f17 + f24 + f26;
	const float kf3 = + f4 + f8 + f9 + f17 + f24 + f26;
	const float kf4 = + f1 + f2 + f8 + f9 + f12 + f15 + f24 + f26;
	const float kf5 = + f6 + f12 + f15 + f17 + f24 + f26;
	const float kf6 = + f4 + f8 + f9 + f17 + f24 + f26;
	const float kf7 = - f8 + f9 - f24 + f26;
	const float kf8 = + f12 - f15 - f24 + f26;
	const float kf9 = + f12 + f15 + f24 + f26;
	const float kf10 = + f8 + f9 + f24 + f26;
	const float kf11 = + f12 - f15 - f24 + f26;
	const float kf12 = - f8 + f9 - f24 + f26;
	const float kf13 = + f8 + f9 + f24 + f26;
	const float kf14 = + f12 + f15 + f24 + f26;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * ux * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f3 = + s0 - s4 - s5 + s14;
	float f5 = + s0 - s4 - s6 + s13;
	float f7 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s11 - 0.5f * s14;
	float f10 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s11 - 0.5f * s14;
	float f11 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s13;
	float f13 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s9 - 0.5f * s14;
	float f14 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s10 - 0.5f * s13;
	float f16 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s13;
	float f18 = - s0 - 0.5f * s2 - 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 + 0.5f * s9 + 0.5f * s10 - 0.5f * s13 - 0.5f * s14;
	float f19 = - 0.25f * s8 + 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f20 = + 0.25f * s7 + 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
	float f21 = - 0.25f * s7 + 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f22 = + 0.25f * s8 + 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
	float f23 = - 0.5f * s1 - 0.5f * s4 - 0.25f * s7 - 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f25 = + 0.5f * s1 - 0.5f * s4 + 0.25f * s7 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
}
else if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f4 + f5 + f6 + f9 + f12 + f14 + f16 + f17 + f20 + f26;
	const float kf1 = + f1 + f9 + f12 + f16 + f20 + f26;
	const float kf2 = - f5 + f6 + f12 - f14 - f16 + f17 - f20 + f26;
	const float kf3 = + f4 + f9 + f14 + f17 + f20 + f26;
	const float kf4 = + f1 + f9 + f12 + f16 + f20 + f26;
	const float kf5 = + f5 + f6 + f12 + f14 + f16 + f17 + f20 + f26;
	const float kf6 = + f4 + f9 + f14 + f17 + f20 + f26;
	const float kf7 = - f14 + f17 - f20 + f26;
	const float kf8 = + f12 - f16 - f20 + f26;
	const float kf9 = + f12 - f16 - f20 + f26;
	const float kf10 = + f12 + f16 + f20 + f26;
	const float kf11 = + f14 + f17 + f20 + f26;
	const float kf12 = - f14 + f17 - f20 + f26;
	const float kf13 = + f14 + f17 + f20 + f26;
	const float kf14 = + f12 + f16 + f20 + f26;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s13;
	float f3 = + s0 - s4 - s5 + s14;
	float f7 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s10 - 0.5f * s14;
	float f8 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s11 - 0.5f * s13;
	float f10 = - s0 - 0.5f * s1 - 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 + 0.5f * s10 + 0.5f * s11 - 0.5f * s13 - 0.5f * s14;
	float f11 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s12 - 0.5f * s13;
	float f13 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s9 - 0.5f * s14;
	float f15 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s12 - 0.5f * s13;
	float f18 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s9 - 0.5f * s14;
	float f19 = - 0.5f * s2 - 0.5f * s5 - 0.25f * s7 - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 - 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f21 = - 0.25f * s7 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13;
	float f22 = + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s14;
	float f23 = - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 + 0.25f * s14;
	float f24 = + 0.25f * s7 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13;
	float f25 = + 0.5f * s2 - 0.5f * s5 + 0.25f * s7 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
}
else if (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f3 + f4 + f6 + f7 + f9 + f12 + f13 + f17 + f22 + f26;
	const float kf1 = + f1 + f7 + f9 + f12 + f22 + f26;
	const float kf2 = + f6 + f12 + f13 + f17 + f22 + f26;
	const float kf3 = - f3 + f4 - f7 + f9 - f13 + f17 - f22 + f26;
	const float kf4 = + f1 + f7 + f9 + f12 + f22 + f26;
	const float kf5 = + f6 + f12 + f13 + f17 + f22 + f26;
	const float kf6 = + f3 + f4 + f7 + f9 + f13 + f17 + f22 + f26;
	const float kf7 = - f13 + f17 - f22 + f26;
	const float kf8 = - f7 + f9 - f22 + f26;
	const float kf9 = - f7 + f9 - f22 + f26;
	const float kf10 = - f13 + f17 - f22 + f26;
	const float kf11 = + f7 + f9 + f22 + f26;
	const float kf12 = + f13 + f17 + f22 + f26;
	const float kf13 = + f13 + f17 + f22 + f26;
	const float kf14 = + f7 + f9 + f22 + f26;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m10 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s13;
	float f5 = + s0 - s4 - s6 + s14;
	float f8 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s10 - 0.5f * s13;
	float f10 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s10 - 0.5f * s13;
	float f11 = - s0 - 0.5f * s1 - 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 + 0.5f * s11 + 0.5f * s12 - 0.5f * s13 - 0.5f * s14;
	float f14 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s9 - 0.5f * s14;
	float f15 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s12 - 0.5f * s13;
	float f16 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s11 - 0.5f * s14;
	float f18 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s9 - 0.5f * s14;
	float f19 = - 0.25f * s7 - 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
	float f20 = + 0.25f * s8 + 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
	float f21 = - 0.5f * s3 - 0.5f * s6 - 0.25f * s7 - 0.25f * s8 + 0.25f * s9 + 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f23 = - 0.25f * s8 - 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
	float f24 = + 0.25f * s7 + 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
	float f25 = + 0.5f * s3 - 0.5f * s6 + 0.25f * s7 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f3 + f5 + f7 + f10 + f11 + f16 + f18 + f23 + f25;
	const float kf1 = + f1 - f2 + f7 - f10 - f11 + f16 + f23 - f25;
	const float kf2 = - f5 - f11 - f16 - f18 - f23 - f25;
	const float kf3 = - f3 - f7 - f10 - f18 - f23 - f25;
	const float kf4 = + f1 + f2 + f7 + f10 + f11 + f16 + f23 + f25;
	const float kf5 = + f5 + f11 + f16 + f18 + f23 + f25;
	const float kf6 = + f3 + f7 + f10 + f18 + f23 + f25;
	const float kf7 = - f7 + f10 - f23 + f25;
	const float kf8 = + f11 - f16 - f23 + f25;
	const float kf9 = - f11 - f16 - f23 - f25;
	const float kf10 = - f7 - f10 - f23 - f25;
	const float kf11 = - f11 + f16 + f23 - f25;
	const float kf12 = + f7 - f10 + f23 - f25;
	const float kf13 = + f7 + f10 + f23 + f25;
	const float kf14 = + f11 + f16 + f23 + f25;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * ux * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f4 = + s0 - s4 - s5 + s14;
	float f6 = + s0 - s4 - s6 + s13;
	float f8 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s11 - 0.5f * s14;
	float f9 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s11 - 0.5f * s14;
	float f12 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s13;
	float f13 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s10 - 0.5f * s13;
	float f14 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s9 - 0.5f * s14;
	float f15 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s13;
	float f17 = - s0 + 0.5f * s2 + 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 - 0.5f * s9 - 0.5f * s10 - 0.5f * s13 - 0.5f * s14;
	float f19 = + 0.25f * s7 - 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f20 = - 0.25f * s8 - 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
	float f21 = + 0.25f * s8 - 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f22 = - 0.25f * s7 - 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
	float f24 = + 0.5f * s1 - 0.5f * s4 - 0.25f * s7 - 0.25f * s8 + 0.25f * s9 + 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f26 = - 0.5f * s1 - 0.5f * s4 + 0.25f * s7 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f3 + f5 + f6 + f10 + f11 + f13 + f15 + f18 + f19 + f25;
	const float kf1 = - f2 - f10 - f11 - f15 - f19 - f25;
	const float kf2 = - f5 + f6 - f11 + f13 + f15 - f18 + f19 - f25;
	const float kf3 = - f3 - f10 - f13 - f18 - f19 - f25;
	const float kf4 = + f2 + f10 + f11 + f15 + f19 + f25;
	const float kf5 = + f5 + f6 + f11 + f13 + f15 + f18 + f19 + f25;
	const float kf6 = + f3 + f10 + f13 + f18 + f19 + f25;
	const float kf7 = - f13 + f18 - f19 + f25;
	const float kf8 = + f11 - f15 - f19 + f25;
	const float kf9 = - f11 + f15 + f19 - f25;
	const float kf10 = - f11 - f15 - f19 - f25;
	const float kf11 = - f13 - f18 - f19 - f25;
	const float kf12 = + f13 - f18 + f19 - f25;
	const float kf13 = + f13 + f18 + f19 + f25;
	const float kf14 = + f11 + f15 + f19 + f25;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s13;
	float f4 = + s0 - s4 - s5 + s14;
	float f7 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s11 - 0.5f * s13;
	float f8 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s10 - 0.5f * s14;
	float f9 = - s0 + 0.5f * s1 + 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 - 0.5f * s10 - 0.5f * s11 - 0.5f * s13 - 0.5f * s14;
	float f12 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s12 - 0.5f * s13;
	float f14 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s9 - 0.5f * s14;
	float f16 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s12 - 0.5f * s13;
	float f17 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s9 - 0.5f * s14;
	float f20 = + 0.5f * s2 - 0.5f * s5 - 0.25f * s7 - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f21 = + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s14;
	float f22 = - 0.25f * s7 - 0.25f * s11 + 0.25f * s12 + 0.25f * s13;
	float f23 = + 0.25f * s7 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13;
	float f24 = - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 + 0.25f * s14;
	float f26 = - 0.5f * s2 - 0.5f * s5 + 0.25f * s7 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
}
else if (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f3 + f4 + f5 + f8 + f10 + f11 + f14 + f18 + f21 + f25;
	const float kf1 = - f2 - f8 - f10 - f11 - f21 - f25;
	const float kf2 = - f5 - f11 - f14 - f18 - f21 - f25;
	const float kf3 = - f3 + f4 + f8 - f10 + f14 - f18 + f21 - f25;
	const float kf4 = + f2 + f8 + f10 + f11 + f21 + f25;
	const float kf5 = + f5 + f11 + f14 + f18 + f21 + f25;
	const float kf6 = + f3 + f4 + f8 + f10 + f14 + f18 + f21 + f25;
	const float kf7 = - f14 + f18 - f21 + f25;
	const float kf8 = - f8 + f10 - f21 + f25;
	const float kf9 = + f8 - f10 + f21 - f25;
	const float kf10 = + f14 - f18 + f21 - f25;
	const float kf11 = - f8 - f10 - f21 - f25;
	const float kf12 = - f14 - f18 - f21 - f25;
	const float kf13 = + f14 + f18 + f21 + f25;
	const float kf14 = + f8 + f10 + f21 + f25;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m10 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s13;
	float f6 = + s0 - s4 - s6 + s14;
	float f7 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s10 - 0.5f * s13;
	float f9 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s10 - 0.5f * s13;
	float f12 = - s0 + 0.5f * s1 + 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 - 0.5f * s11 - 0.5f * s12 - 0.5f * s13 - 0.5f * s14;
	float f13 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s9 - 0.5f * s14;
	float f15 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s11 - 0.5f * s14;
	float f16 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s12 - 0.5f * s13;
	float f17 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s9 - 0.5f * s14;
	float f19 = + 0.25f * s8 - 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f20 = - 0.25f * s7 + 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f22 = + 0.5f * s3 - 0.5f * s6 - 0.25f * s7 - 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f23 = + 0.25f * s7 - 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f24 = - 0.25f * s8 + 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f26 = - 0.5f * s3 - 0.5f * s6 + 0.25f * s7 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
}
else if (outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f3 + f6 + f7 + f10 + f12 + f13 + f15 + f19 + f22;
	const float kf1 = + f1 - f2 + f7 - f10 + f12 - f15 - f19 + f22;
	const float kf2 = + f6 + f12 + f13 + f15 + f19 + f22;
	const float kf3 = - f3 - f7 - f10 - f13 - f19 - f22;
	const float kf4 = + f1 + f2 + f7 + f10 + f12 + f15 + f19 + f22;
	const float kf5 = + f6 + f12 + f13 + f15 + f19 + f22;
	const float kf6 = + f3 + f7 + f10 + f13 + f19 + f22;
	const float kf7 = - f7 + f10 + f19 - f22;
	const float kf8 = + f12 - f15 - f19 + f22;
	const float kf9 = + f12 + f15 + f19 + f22;
	const float kf10 = - f7 - f10 - f19 - f22;
	const float kf11 = + f12 - f15 - f19 + f22;
	const float kf12 = + f7 - f10 - f19 + f22;
	const float kf13 = + f7 + f10 + f19 + f22;
	const float kf14 = + f12 + f15 + f19 + f22;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * ux * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f4 = + s0 - s4 - s5 + s14;
	float f5 = + s0 - s4 - s6 + s13;
	float f8 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s11 - 0.5f * s14;
	float f9 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s11 - 0.5f * s14;
	float f11 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s13;
	float f14 = - s0 - 0.5f * s2 + 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 + 0.5f * s9 - 0.5f * s10 - 0.5f * s13 - 0.5f * s14;
	float f16 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s13;
	float f17 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s9 - 0.5f * s14;
	float f18 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s10 - 0.5f * s13;
	float f20 = - 0.5f * s1 - 0.5f * s4 + 0.25f * s7 - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f21 = + 0.5f * s1 - 0.5f * s4 - 0.25f * s7 + 0.25f * s8 - 0.25f * s9 + 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f23 = - 0.25f * s7 - 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
	float f24 = - 0.25f * s8 + 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f25 = + 0.25f * s7 - 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f26 = + 0.25f * s8 + 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
}
else if (outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f3 + f5 + f6 + f7 + f12 + f13 + f16 + f18 + f22 + f23;
	const float kf1 = + f1 + f7 + f12 + f16 + f22 + f23;
	const float kf2 = - f5 + f6 + f12 + f13 - f16 - f18 + f22 - f23;
	const float kf3 = - f3 - f7 - f13 - f18 - f22 - f23;
	const float kf4 = + f1 + f7 + f12 + f16 + f22 + f23;
	const float kf5 = + f5 + f6 + f12 + f13 + f16 + f18 + f22 + f23;
	const float kf6 = + f3 + f7 + f13 + f18 + f22 + f23;
	const float kf7 = - f13 + f18 - f22 + f23;
	const float kf8 = + f12 - f16 + f22 - f23;
	const float kf9 = + f12 - f16 + f22 - f23;
	const float kf10 = + f12 + f16 + f22 + f23;
	const float kf11 = - f13 - f18 - f22 - f23;
	const float kf12 = + f13 - f18 + f22 - f23;
	const float kf13 = + f13 + f18 + f22 + f23;
	const float kf14 = + f12 + f16 + f22 + f23;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s13;
	float f4 = + s0 - s4 - s5 + s14;
	float f8 = - s0 - 0.5f * s1 + 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 + 0.5f * s10 - 0.5f * s11 - 0.5f * s13 - 0.5f * s14;
	float f9 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s10 - 0.5f * s14;
	float f10 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s11 - 0.5f * s13;
	float f11 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s12 - 0.5f * s13;
	float f14 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s9 - 0.5f * s14;
	float f15 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s12 - 0.5f * s13;
	float f17 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s9 - 0.5f * s14;
	float f19 = - 0.25f * s7 - 0.25f * s11 + 0.25f * s12 + 0.25f * s13;
	float f20 = - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 + 0.25f * s14;
	float f21 = + 0.5f * s2 - 0.5f * s5 - 0.25f * s7 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f24 = - 0.5f * s2 - 0.5f * s5 + 0.25f * s7 - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f25 = + 0.25f * s7 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13;
	float f26 = + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s14;
}
else if (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 0)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f3 + f4 + f5 + f7 + f9 + f14 + f16 + f18 + f20 + f23;
	const float kf1 = + f1 + f7 + f9 + f16 + f20 + f23;
	const float kf2 = - f5 - f14 - f16 - f18 - f20 - f23;
	const float kf3 = - f3 + f4 - f7 + f9 + f14 - f18 + f20 - f23;
	const float kf4 = + f1 + f7 + f9 + f16 + f20 + f23;
	const float kf5 = + f5 + f14 + f16 + f18 + f20 + f23;
	const float kf6 = + f3 + f4 + f7 + f9 + f14 + f18 + f20 + f23;
	const float kf7 = - f14 + f18 - f20 + f23;
	const float kf8 = - f7 + f9 + f20 - f23;
	const float kf9 = - f7 + f9 + f20 - f23;
	const float kf10 = + f14 - f18 + f20 - f23;
	const float kf11 = + f7 + f9 + f20 + f23;
	const float kf12 = - f14 - f18 - f20 - f23;
	const float kf13 = + f14 + f18 + f20 + f23;
	const float kf14 = + f7 + f9 + f20 + f23;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m10 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s13;
	float f6 = + s0 - s4 - s6 + s14;
	float f8 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s10 - 0.5f * s13;
	float f10 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s10 - 0.5f * s13;
	float f11 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s12 - 0.5f * s13;
	float f12 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s11 - 0.5f * s14;
	float f13 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s9 - 0.5f * s14;
	float f15 = - s0 - 0.5f * s1 + 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 + 0.5f * s11 - 0.5f * s12 - 0.5f * s13 - 0.5f * s14;
	float f17 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s9 - 0.5f * s14;
	float f19 = + 0.5f * s3 - 0.5f * s6 - 0.25f * s7 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 - 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f21 = - 0.25f * s7 + 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f22 = - 0.25f * s8 - 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
	float f24 = - 0.5f * s3 - 0.5f * s6 + 0.25f * s7 - 0.25f * s8 + 0.25f * s9 + 0.25f * s10 - 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f25 = + 0.25f * s7 - 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f26 = + 0.25f * s8 + 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
}
else if (outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f2 + f4 + f5 + f8 + f9 + f11 + f14 + f16 + f20 + f21;
	const float kf1 = + f1 - f2 - f8 + f9 - f11 + f16 + f20 - f21;
	const float kf2 = - f5 - f11 - f14 - f16 - f20 - f21;
	const float kf3 = + f4 + f8 + f9 + f14 + f20 + f21;
	const float kf4 = + f1 + f2 + f8 + f9 + f11 + f16 + f20 + f21;
	const float kf5 = + f5 + f11 + f14 + f16 + f20 + f21;
	const float kf6 = + f4 + f8 + f9 + f14 + f20 + f21;
	const float kf7 = - f8 + f9 + f20 - f21;
	const float kf8 = + f11 - f16 - f20 + f21;
	const float kf9 = - f11 - f16 - f20 - f21;
	const float kf10 = + f8 + f9 + f20 + f21;
	const float kf11 = - f11 + f16 + f20 - f21;
	const float kf12 = - f8 + f9 + f20 - f21;
	const float kf13 = + f8 + f9 + f20 + f21;
	const float kf14 = + f11 + f16 + f20 + f21;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * ux * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f3 = + s0 - s4 - s5 + s14;
	float f6 = + s0 - s4 - s6 + s13;
	float f7 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s11 - 0.5f * s14;
	float f10 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s11 - 0.5f * s14;
	float f12 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s13;
	float f13 = - s0 + 0.5f * s2 - 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 - 0.5f * s9 + 0.5f * s10 - 0.5f * s13 - 0.5f * s14;
	float f15 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s13;
	float f17 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s10 - 0.5f * s13;
	float f18 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s9 - 0.5f * s14;
	float f19 = + 0.5f * s1 - 0.5f * s4 + 0.25f * s7 - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f22 = - 0.5f * s1 - 0.5f * s4 - 0.25f * s7 + 0.25f * s8 + 0.25f * s9 - 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f23 = - 0.25f * s8 - 0.25f * s9 + 0.25f * s11 + 0.25f * s14;
	float f24 = - 0.25f * s7 + 0.25f * s10 - 0.25f * s12 + 0.25f * s13;
	float f25 = + 0.25f * s8 - 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f26 = + 0.25f * s7 + 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
}
else if (outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f4 + f5 + f6 + f8 + f11 + f14 + f15 + f17 + f21 + f24;
	const float kf1 = - f2 - f8 - f11 - f15 - f21 - f24;
	const float kf2 = - f5 + f6 - f11 - f14 + f15 + f17 - f21 + f24;
	const float kf3 = + f4 + f8 + f14 + f17 + f21 + f24;
	const float kf4 = + f2 + f8 + f11 + f15 + f21 + f24;
	const float kf5 = + f5 + f6 + f11 + f14 + f15 + f17 + f21 + f24;
	const float kf6 = + f4 + f8 + f14 + f17 + f21 + f24;
	const float kf7 = - f14 + f17 - f21 + f24;
	const float kf8 = + f11 - f15 + f21 - f24;
	const float kf9 = - f11 + f15 - f21 + f24;
	const float kf10 = - f11 - f15 - f21 - f24;
	const float kf11 = + f14 + f17 + f21 + f24;
	const float kf12 = - f14 + f17 - f21 + f24;
	const float kf13 = + f14 + f17 + f21 + f24;
	const float kf14 = + f11 + f15 + f21 + f24;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uy;
	const float m9 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m10 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s13;
	float f3 = + s0 - s4 - s5 + s14;
	float f7 = - s0 + 0.5f * s1 - 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 - 0.5f * s10 + 0.5f * s11 - 0.5f * s13 - 0.5f * s14;
	float f9 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s11 - 0.5f * s13;
	float f10 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s10 - 0.5f * s14;
	float f12 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s12 - 0.5f * s13;
	float f13 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s9 - 0.5f * s14;
	float f16 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s12 - 0.5f * s13;
	float f18 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s9 - 0.5f * s14;
	float f19 = - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 + 0.25f * s14;
	float f20 = - 0.25f * s7 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13;
	float f22 = - 0.5f * s2 - 0.5f * s5 - 0.25f * s7 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 - 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f23 = + 0.5f * s2 - 0.5f * s5 + 0.25f * s7 - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 - 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f25 = + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s14;
	float f26 = + 0.25f * s7 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13;
}
else if (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 0)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f3 + f4 + f6 + f8 + f10 + f13 + f15 + f17 + f19 + f24;
	const float kf1 = - f2 - f8 - f10 - f15 - f19 - f24;
	const float kf2 = + f6 + f13 + f15 + f17 + f19 + f24;
	const float kf3 = - f3 + f4 + f8 - f10 - f13 + f17 - f19 + f24;
	const float kf4 = + f2 + f8 + f10 + f15 + f19 + f24;
	const float kf5 = + f6 + f13 + f15 + f17 + f19 + f24;
	const float kf6 = + f3 + f4 + f8 + f10 + f13 + f17 + f19 + f24;
	const float kf7 = - f13 + f17 - f19 + f24;
	const float kf8 = - f8 + f10 + f19 - f24;
	const float kf9 = + f8 - f10 - f19 + f24;
	const float kf10 = - f13 + f17 - f19 + f24;
	const float kf11 = - f8 - f10 - f19 - f24;
	const float kf12 = + f13 + f17 + f19 + f24;
	const float kf13 = + f13 + f17 + f19 + f24;
	const float kf14 = + f8 + f10 + f19 + f24;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m10 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m11 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m12 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m13 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m14 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s13;
	float f5 = + s0 - s4 - s6 + s14;
	float f7 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s10 - 0.5f * s13;
	float f9 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s10 - 0.5f * s13;
	float f11 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s11 - 0.5f * s14;
	float f12 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s12 - 0.5f * s13;
	float f14 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s9 - 0.5f * s14;
	float f16 = - s0 + 0.5f * s1 - 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 - 0.5f * s11 + 0.5f * s12 - 0.5f * s13 - 0.5f * s14;
	float f18 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s9 - 0.5f * s14;
	float f20 = - 0.5f * s3 - 0.5f * s6 - 0.25f * s7 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f21 = - 0.25f * s8 + 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f22 = - 0.25f * s7 - 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
	float f23 = + 0.5f * s3 - 0.5f * s6 + 0.25f * s7 - 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13 + 0.25f * s14;
	float f25 = + 0.25f * s8 - 0.25f * s9 - 0.25f * s11 + 0.25f * s14;
	float f26 = + 0.25f * s7 + 0.25f * s10 + 0.25f * s12 + 0.25f * s13;
}
else if (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f4 + f6 + f9 + f12 + f17 + f26;
	const float kf1 = + f1 + f9 + f12 + f26;
	const float kf2 = + f6 + f12 + f17 + f26;
	const float kf3 = + f4 + f9 + f17 + f26;
	const float kf4 = + f1 + f9 + f12 + f26;
	const float kf5 = + f6 + f12 + f17 + f26;
	const float kf6 = + f4 + f9 + f17 + f26;
	const float kf7 = + f17 + f26;
	const float kf8 = + f9 + f26;
	const float kf9 = + f12 + f26;
	const float kf10 = + f12 + f26;
	const float kf11 = + f9 + f26;
	const float kf12 = + f12 + f26;
	const float kf13 = + f17 + f26;
	const float kf14 = + f9 + f26;
	const float kf15 = + f17 + f26;
	const float kf16 = + f17 + f26;
	const float kf17 = + f9 + f26;
	const float kf18 = + f12 + f26;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s16;
	float f3 = + s0 - s4 - s5 + s18;
	float f5 = + s0 - s4 - s6 + s17;
	float f7 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s18;
	float f8 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s13 - 0.5f * s16;
	float f10 = - s0 - 0.5f * s1 - 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 + 0.5f * s12 + 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f11 = - s0 - 0.5f * s1 - 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 + 0.5f * s14 + 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f13 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s10 - 0.5f * s18;
	float f14 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s11 - 0.5f * s17;
	float f15 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s15 - 0.5f * s16;
	float f16 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s14 - 0.5f * s17;
	float f18 = - s0 - 0.5f * s2 - 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 + 0.5f * s10 + 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f19 = - 0.5f * s2 - 0.5f * s5 - 0.25f * s7 - 0.25f * s9 + 0.25f * s10 - 0.25f * s12 - 0.25f * s13 + 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f20 = + 0.25f * s8 + 0.25f * s11 + 0.25f * s14 + 0.25f * s17;
	float f21 = - 0.5f * s3 - 0.5f * s6 - 0.25f * s7 - 0.25f * s8 + 0.25f * s11 + 0.25f * s13 - 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f22 = + 0.25f * s9 + 0.25f * s10 + 0.25f * s12 + 0.25f * s18;
	float f23 = - 0.5f * s1 - 0.5f * s4 - 0.25f * s8 - 0.25f * s9 - 0.25f * s10 - 0.25f * s11 + 0.25f * s12 + 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f24 = + 0.25f * s7 + 0.25f * s13 + 0.25f * s15 + 0.25f * s16;
	float f25 = + s0 + 0.5f * s1 + 0.5f * s2 + 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 + 0.25f * s7 + 0.25f * s8 + 0.25f * s9 - 0.25f * s10 - 0.25f * s11 - 0.25f * s12 - 0.25f * s13 - 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
}
else if (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f3 + f5 + f10 + f11 + f18 + f25;
	const float kf1 = - f2 - f10 - f11 - f25;
	const float kf2 = - f5 - f11 - f18 - f25;
	const float kf3 = - f3 - f10 - f18 - f25;
	const float kf4 = + f2 + f10 + f11 + f25;
	const float kf5 = + f5 + f11 + f18 + f25;
	const float kf6 = + f3 + f10 + f18 + f25;
	const float kf7 = + f18 + f25;
	const float kf8 = + f10 + f25;
	const float kf9 = + f11 + f25;
	const float kf10 = - f11 - f25;
	const float kf11 = - f10 - f25;
	const float kf12 = - f11 - f25;
	const float kf13 = - f18 - f25;
	const float kf14 = - f10 - f25;
	const float kf15 = - f18 - f25;
	const float kf16 = + f18 + f25;
	const float kf17 = + f10 + f25;
	const float kf18 = + f11 + f25;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s16;
	float f4 = + s0 - s4 - s5 + s18;
	float f6 = + s0 - s4 - s6 + s17;
	float f7 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s13 - 0.5f * s16;
	float f8 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s18;
	float f9 = - s0 + 0.5f * s1 + 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 - 0.5f * s12 - 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f12 = - s0 + 0.5f * s1 + 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 - 0.5f * s14 - 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f13 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s11 - 0.5f * s17;
	float f14 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s10 - 0.5f * s18;
	float f15 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s14 - 0.5f * s17;
	float f16 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s15 - 0.5f * s16;
	float f17 = - s0 + 0.5f * s2 + 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 - 0.5f * s10 - 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f19 = + 0.25f * s8 - 0.25f * s11 - 0.25f * s14 + 0.25f * s17;
	float f20 = + 0.5f * s2 - 0.5f * s5 - 0.25f * s7 - 0.25f * s9 - 0.25f * s10 + 0.25f * s12 + 0.25f * s13 - 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f21 = + 0.25f * s9 - 0.25f * s10 - 0.25f * s12 + 0.25f * s18;
	float f22 = + 0.5f * s3 - 0.5f * s6 - 0.25f * s7 - 0.25f * s8 - 0.25f * s11 - 0.25f * s13 + 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f23 = + 0.25f * s7 - 0.25f * s13 - 0.25f * s15 + 0.25f * s16;
	float f24 = + 0.5f * s1 - 0.5f * s4 - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 + 0.25f * s11 - 0.25f * s12 - 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f26 = + s0 - 0.5f * s1 - 0.5f * s2 - 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 + 0.25f * s7 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
}
else if (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f4 + f6 + f8 + f15 + f17 + f24;
	const float kf1 = - f2 - f8 - f15 - f24;
	const float kf2 = + f6 + f15 + f17 + f24;
	const float kf3 = + f4 + f8 + f17 + f24;
	const float kf4 = + f2 + f8 + f15 + f24;
	const float kf5 = + f6 + f15 + f17 + f24;
	const float kf6 = + f4 + f8 + f17 + f24;
	const float kf7 = + f17 + f24;
	const float kf8 = - f8 - f24;
	const float kf9 = - f15 - f24;
	const float kf10 = + f15 + f24;
	const float kf11 = + f8 + f24;
	const float kf12 = - f15 - f24;
	const float kf13 = + f17 + f24;
	const float kf14 = - f8 - f24;
	const float kf15 = + f17 + f24;
	const float kf16 = + f17 + f24;
	const float kf17 = + f8 + f24;
	const float kf18 = + f15 + f24;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s16;
	float f3 = + s0 - s4 - s5 + s18;
	float f5 = + s0 - s4 - s6 + s17;
	float f7 = - s0 + 0.5f * s1 - 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 - 0.5f * s12 + 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f9 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s13 - 0.5f * s16;
	float f10 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s18;
	float f11 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s14 - 0.5f * s17;
	float f12 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s15 - 0.5f * s16;
	float f13 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s10 - 0.5f * s18;
	float f14 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s11 - 0.5f * s17;
	float f16 = - s0 + 0.5f * s1 - 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 - 0.5f * s14 + 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f18 = - s0 - 0.5f * s2 - 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 + 0.5f * s10 + 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f19 = - 0.25f * s9 + 0.25f * s10 - 0.25f * s12 + 0.25f * s18;
	float f20 = - 0.5f * s3 - 0.5f * s6 - 0.25f * s7 + 0.25f * s8 + 0.25f * s11 + 0.25f * s13 + 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f21 = - 0.25f * s8 + 0.25f * s11 - 0.25f * s14 + 0.25f * s17;
	float f22 = - 0.5f * s2 - 0.5f * s5 - 0.25f * s7 + 0.25f * s9 + 0.25f * s10 + 0.25f * s12 - 0.25f * s13 + 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f23 = + s0 - 0.5f * s1 + 0.5f * s2 + 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 + 0.25f * s7 - 0.25f * s8 - 0.25f * s9 - 0.25f * s10 - 0.25f * s11 + 0.25f * s12 - 0.25f * s13 + 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
	float f25 = + 0.5f * s1 - 0.5f * s4 + 0.25f * s8 + 0.25f * s9 - 0.25f * s10 - 0.25f * s11 - 0.25f * s12 - 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f26 = + 0.25f * s7 + 0.25f * s13 + 0.25f * s15 + 0.25f * s16;
}
else if (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f4 + f5 + f9 + f14 + f16 + f20;
	const float kf1 = + f1 + f9 + f16 + f20;
	const float kf2 = - f5 - f14 - f16 - f20;
	const float kf3 = + f4 + f9 + f14 + f20;
	const float kf4 = + f1 + f9 + f16 + f20;
	const float kf5 = + f5 + f14 + f16 + f20;
	const float kf6 = + f4 + f9 + f14 + f20;
	const float kf7 = - f14 - f20;
	const float kf8 = + f9 + f20;
	const float kf9 = - f16 - f20;
	const float kf10 = - f16 - f20;
	const float kf11 = + f9 + f20;
	const float kf12 = + f16 + f20;
	const float kf13 = + f14 + f20;
	const float kf14 = + f9 + f20;
	const float kf15 = - f14 - f20;
	const float kf16 = + f14 + f20;
	const float kf17 = + f9 + f20;
	const float kf18 = + f16 + f20;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s16;
	float f3 = + s0 - s4 - s5 + s18;
	float f6 = + s0 - s4 - s6 + s17;
	float f7 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s18;
	float f8 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s13 - 0.5f * s16;
	float f10 = - s0 - 0.5f * s1 - 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 + 0.5f * s12 + 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f11 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s15 - 0.5f * s16;
	float f12 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s14 - 0.5f * s17;
	float f13 = - s0 + 0.5f * s2 - 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 - 0.5f * s10 + 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f15 = - s0 - 0.5f * s1 + 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 + 0.5f * s14 - 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f17 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s11 - 0.5f * s17;
	float f18 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s10 - 0.5f * s18;
	float f19 = + s0 + 0.5f * s1 - 0.5f * s2 + 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 - 0.25f * s7 + 0.25f * s8 - 0.25f * s9 + 0.25f * s10 - 0.25f * s11 - 0.25f * s12 - 0.25f * s13 - 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
	float f21 = - 0.25f * s7 + 0.25f * s13 - 0.25f * s15 + 0.25f * s16;
	float f22 = - 0.5f * s1 - 0.5f * s4 - 0.25f * s8 + 0.25f * s9 + 0.25f * s10 - 0.25f * s11 + 0.25f * s12 + 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f23 = - 0.25f * s9 - 0.25f * s10 + 0.25f * s12 + 0.25f * s18;
	float f24 = - 0.5f * s3 - 0.5f * s6 + 0.25f * s7 - 0.25f * s8 + 0.25f * s11 + 0.25f * s13 - 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f25 = + 0.5f * s2 - 0.5f * s5 + 0.25f * s7 + 0.25f * s9 - 0.25f * s10 - 0.25f * s12 - 0.25f * s13 - 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f26 = + 0.25f * s8 + 0.25f * s11 + 0.25f * s14 + 0.25f * s17;
}
else if (outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f3 + f6 + f7 + f12 + f13 + f22;
	const float kf1 = + f1 + f7 + f12 + f22;
	const float kf2 = + f6 + f12 + f13 + f22;
	const float kf3 = - f3 - f7 - f13 - f22;
	const float kf4 = + f1 + f7 + f12 + f22;
	const float kf5 = + f6 + f12 + f13 + f22;
	const float kf6 = + f3 + f7 + f13 + f22;
	const float kf7 = - f13 - f22;
	const float kf8 = - f7 - f22;
	const float kf9 = + f12 + f22;
	const float kf10 = + f12 + f22;
	const float kf11 = - f7 - f22;
	const float kf12 = + f12 + f22;
	const float kf13 = - f13 - f22;
	const float kf14 = + f7 + f22;
	const float kf15 = + f13 + f22;
	const float kf16 = + f13 + f22;
	const float kf17 = + f7 + f22;
	const float kf18 = + f12 + f22;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s16;
	float f4 = + s0 - s4 - s5 + s18;
	float f5 = + s0 - s4 - s6 + s17;
	float f8 = - s0 - 0.5f * s1 + 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 + 0.5f * s12 - 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f9 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s18;
	float f10 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s13 - 0.5f * s16;
	float f11 = - s0 - 0.5f * s1 - 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 + 0.5f * s14 + 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f14 = - s0 - 0.5f * s2 + 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 + 0.5f * s10 - 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f15 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s15 - 0.5f * s16;
	float f16 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s14 - 0.5f * s17;
	float f17 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s10 - 0.5f * s18;
	float f18 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s11 - 0.5f * s17;
	float f19 = - 0.25f * s7 - 0.25f * s13 + 0.25f * s15 + 0.25f * s16;
	float f20 = - 0.5f * s1 - 0.5f * s4 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f21 = + s0 + 0.5f * s1 + 0.5f * s2 - 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 - 0.25f * s7 - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13 - 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
	float f23 = - 0.25f * s8 - 0.25f * s11 + 0.25f * s14 + 0.25f * s17;
	float f24 = - 0.5f * s2 - 0.5f * s5 + 0.25f * s7 - 0.25f * s9 + 0.25f * s10 - 0.25f * s12 + 0.25f * s13 + 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f25 = + 0.5f * s3 - 0.5f * s6 + 0.25f * s7 + 0.25f * s8 - 0.25f * s11 - 0.25f * s13 - 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f26 = + 0.25f * s9 + 0.25f * s10 + 0.25f * s12 + 0.25f * s18;
}
else if (outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f1 + f3 + f5 + f7 + f16 + f18 + f23;
	const float kf1 = + f1 + f7 + f16 + f23;
	const float kf2 = - f5 - f16 - f18 - f23;
	const float kf3 = - f3 - f7 - f18 - f23;
	const float kf4 = + f1 + f7 + f16 + f23;
	const float kf5 = + f5 + f16 + f18 + f23;
	const float kf6 = + f3 + f7 + f18 + f23;
	const float kf7 = + f18 + f23;
	const float kf8 = - f7 - f23;
	const float kf9 = - f16 - f23;
	const float kf10 = - f16 - f23;
	const float kf11 = - f7 - f23;
	const float kf12 = + f16 + f23;
	const float kf13 = - f18 - f23;
	const float kf14 = + f7 + f23;
	const float kf15 = - f18 - f23;
	const float kf16 = + f18 + f23;
	const float kf17 = + f7 + f23;
	const float kf18 = + f16 + f23;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f2 = + s0 - s5 - s6 + s16;
	float f4 = + s0 - s4 - s5 + s18;
	float f6 = + s0 - s4 - s6 + s17;
	float f8 = - s0 - 0.5f * s1 + 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 + 0.5f * s12 - 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f9 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s12 - 0.5f * s18;
	float f10 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s13 - 0.5f * s16;
	float f11 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s15 - 0.5f * s16;
	float f12 = + 0.5f * s1 + 0.5f * s4 - 0.5f * s14 - 0.5f * s17;
	float f13 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s11 - 0.5f * s17;
	float f14 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s10 - 0.5f * s18;
	float f15 = - s0 - 0.5f * s1 + 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 + 0.5f * s14 - 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f17 = - s0 + 0.5f * s2 + 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 - 0.5f * s10 - 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f19 = + 0.5f * s3 - 0.5f * s6 - 0.25f * s7 + 0.25f * s8 - 0.25f * s11 - 0.25f * s13 - 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f20 = - 0.25f * s9 - 0.25f * s10 + 0.25f * s12 + 0.25f * s18;
	float f21 = + 0.5f * s2 - 0.5f * s5 - 0.25f * s7 + 0.25f * s9 - 0.25f * s10 - 0.25f * s12 + 0.25f * s13 - 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f22 = - 0.25f * s8 - 0.25f * s11 + 0.25f * s14 + 0.25f * s17;
	float f24 = + s0 + 0.5f * s1 - 0.5f * s2 - 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 + 0.25f * s7 - 0.25f * s8 - 0.25f * s9 + 0.25f * s10 + 0.25f * s11 - 0.25f * s12 + 0.25f * s13 - 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
	float f25 = + 0.25f * s7 - 0.25f * s13 - 0.25f * s15 + 0.25f * s16;
	float f26 = - 0.5f * s1 - 0.5f * s4 + 0.25f * s8 + 0.25f * s9 + 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
}
else if (outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == -1)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f3 + f6 + f10 + f13 + f15 + f19;
	const float kf1 = - f2 - f10 - f15 - f19;
	const float kf2 = + f6 + f13 + f15 + f19;
	const float kf3 = - f3 - f10 - f13 - f19;
	const float kf4 = + f2 + f10 + f15 + f19;
	const float kf5 = + f6 + f13 + f15 + f19;
	const float kf6 = + f3 + f10 + f13 + f19;
	const float kf7 = - f13 - f19;
	const float kf8 = + f10 + f19;
	const float kf9 = - f15 - f19;
	const float kf10 = + f15 + f19;
	const float kf11 = - f10 - f19;
	const float kf12 = - f15 - f19;
	const float kf13 = - f13 - f19;
	const float kf14 = - f10 - f19;
	const float kf15 = + f13 + f19;
	const float kf16 = + f13 + f19;
	const float kf17 = + f10 + f19;
	const float kf18 = + f15 + f19;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s16;
	float f4 = + s0 - s4 - s5 + s18;
	float f5 = + s0 - s4 - s6 + s17;
	float f7 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s13 - 0.5f * s16;
	float f8 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s18;
	float f9 = - s0 + 0.5f * s1 + 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 - 0.5f * s12 - 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f11 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s14 - 0.5f * s17;
	float f12 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s15 - 0.5f * s16;
	float f14 = - s0 - 0.5f * s2 + 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 + 0.5f * s10 - 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f16 = - s0 + 0.5f * s1 - 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 - 0.5f * s14 + 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f17 = + 0.5f * s2 + 0.5f * s5 - 0.5f * s10 - 0.5f * s18;
	float f18 = - 0.5f * s3 + 0.5f * s6 + 0.5f * s11 - 0.5f * s17;
	float f20 = + s0 - 0.5f * s1 + 0.5f * s2 - 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 - 0.25f * s7 + 0.25f * s8 - 0.25f * s9 - 0.25f * s10 + 0.25f * s11 + 0.25f * s12 + 0.25f * s13 + 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
	float f21 = + 0.5f * s1 - 0.5f * s4 - 0.25f * s8 + 0.25f * s9 - 0.25f * s10 + 0.25f * s11 - 0.25f * s12 - 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f22 = - 0.25f * s7 - 0.25f * s13 + 0.25f * s15 + 0.25f * s16;
	float f23 = + 0.5f * s3 - 0.5f * s6 + 0.25f * s7 - 0.25f * s8 - 0.25f * s11 - 0.25f * s13 + 0.25f * s14 - 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
	float f24 = - 0.25f * s9 + 0.25f * s10 - 0.25f * s12 + 0.25f * s18;
	float f25 = + 0.25f * s8 - 0.25f * s11 - 0.25f * s14 + 0.25f * s17;
	float f26 = - 0.5f * s2 - 0.5f * s5 + 0.25f * s7 + 0.25f * s9 + 0.25f * s10 + 0.25f * s12 + 0.25f * s13 + 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
}
else if (outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 1)
{
	// Multiply K fk
	const float kf0 = + f0 + f2 + f4 + f5 + f8 + f11 + f14 + f21;
	const float kf1 = - f2 - f8 - f11 - f21;
	const float kf2 = - f5 - f11 - f14 - f21;
	const float kf3 = + f4 + f8 + f14 + f21;
	const float kf4 = + f2 + f8 + f11 + f21;
	const float kf5 = + f5 + f11 + f14 + f21;
	const float kf6 = + f4 + f8 + f14 + f21;
	const float kf7 = - f14 - f21;
	const float kf8 = - f8 - f21;
	const float kf9 = + f11 + f21;
	const float kf10 = - f11 - f21;
	const float kf11 = + f8 + f21;
	const float kf12 = - f11 - f21;
	const float kf13 = + f14 + f21;
	const float kf14 = - f8 - f21;
	const float kf15 = - f14 - f21;
	const float kf16 = + f14 + f21;
	const float kf17 = + f8 + f21;
	const float kf18 = + f11 + f21;

	// Calculate equilibrium moments
	const float m0 = rho;
	const float m1 = rho * ux;
	const float m2 = rho * uy;
	const float m3 = rho * uz;
	const float m4 = (1.f/3.f) * rho + rho * ux * ux;
	const float m5 = (1.f/3.f) * rho + rho * uy * uy;
	const float m6 = (1.f/3.f) * rho + rho * uz * uz;
	const float m7 = rho * uy * uz;
	const float m8 = rho * ux * uz;
	const float m9 = rho * ux * uy;
	const float m10 = (1.f/3.f) * rho * uy + rho * ux * ux * uy;
	const float m11 = (1.f/3.f) * rho * uz + rho * ux * ux * uz;
	const float m12 = (1.f/3.f) * rho * ux + rho * ux * uy * uy;
	const float m13 = (1.f/3.f) * rho * uz + rho * uy * uy * uz;
	const float m14 = (1.f/3.f) * rho * ux + rho * ux * uz * uz;
	const float m15 = (1.f/3.f) * rho * uy + rho * uy * uz * uz;
	const float m16 = (1.f/9.f) * rho + (1.f/3.f) * rho * uy * uy + (1.f/3.f) * rho * uz * uz + rho * uy * uy * uz * uz;
	const float m17 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uz * uz + rho * ux * ux * uz * uz;
	const float m18 = (1.f/9.f) * rho + (1.f/3.f) * rho * ux * ux + (1.f/3.f) * rho * uy * uy + rho * ux * ux * uy * uy;

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
	const float s9 = m9 - kf9;
	const float s10 = m10 - kf10;
	const float s11 = m11 - kf11;
	const float s12 = m12 - kf12;
	const float s13 = m13 - kf13;
	const float s14 = m14 - kf14;
	const float s15 = m15 - kf15;
	const float s16 = m16 - kf16;
	const float s17 = m17 - kf17;
	const float s18 = m18 - kf18;

	// Multiply U^-1 * (m - Kfk) to get unknown distributions
	float f1 = + s0 - s5 - s6 + s16;
	float f3 = + s0 - s4 - s5 + s18;
	float f6 = + s0 - s4 - s6 + s17;
	float f7 = - s0 + 0.5f * s1 - 0.5f * s3 + 0.5f * s4 + s5 + 0.5f * s6 - 0.5f * s12 + 0.5f * s13 - 0.5f * s16 - 0.5f * s18;
	float f9 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s13 - 0.5f * s16;
	float f10 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s12 - 0.5f * s18;
	float f12 = - s0 + 0.5f * s1 + 0.5f * s2 + 0.5f * s4 + 0.5f * s5 + s6 - 0.5f * s14 - 0.5f * s15 - 0.5f * s16 - 0.5f * s17;
	float f13 = - s0 + 0.5f * s2 - 0.5f * s3 + s4 + 0.5f * s5 + 0.5f * s6 - 0.5f * s10 + 0.5f * s11 - 0.5f * s17 - 0.5f * s18;
	float f15 = - 0.5f * s1 + 0.5f * s4 + 0.5f * s14 - 0.5f * s17;
	float f16 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s15 - 0.5f * s16;
	float f17 = + 0.5f * s3 + 0.5f * s6 - 0.5f * s11 - 0.5f * s17;
	float f18 = - 0.5f * s2 + 0.5f * s5 + 0.5f * s10 - 0.5f * s18;
	float f19 = + 0.5f * s1 - 0.5f * s4 + 0.25f * s8 - 0.25f * s9 + 0.25f * s10 - 0.25f * s11 - 0.25f * s12 - 0.25f * s14 + 0.25f * s17 + 0.25f * s18;
	float f20 = - 0.25f * s7 + 0.25f * s13 - 0.25f * s15 + 0.25f * s16;
	float f22 = + s0 - 0.5f * s1 - 0.5f * s2 + 0.5f * s3 - 0.5f * s4 - 0.5f * s5 - 0.5f * s6 - 0.25f * s7 - 0.25f * s8 + 0.25f * s9 + 0.25f * s10 - 0.25f * s11 + 0.25f * s12 - 0.25f * s13 + 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17 + 0.25f * s18;
	float f23 = + 0.5f * s2 - 0.5f * s5 + 0.25f * s7 - 0.25f * s9 - 0.25f * s10 + 0.25f * s12 - 0.25f * s13 - 0.25f * s15 + 0.25f * s16 + 0.25f * s18;
	float f24 = - 0.25f * s8 + 0.25f * s11 - 0.25f * s14 + 0.25f * s17;
	float f25 = + 0.25f * s9 - 0.25f * s10 - 0.25f * s12 + 0.25f * s18;
	float f26 = - 0.5f * s3 - 0.5f * s6 + 0.25f * s7 + 0.25f * s8 + 0.25f * s11 + 0.25f * s13 + 0.25f * s14 + 0.25f * s15 + 0.25f * s16 + 0.25f * s17;
}