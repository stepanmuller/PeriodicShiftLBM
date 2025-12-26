__host__ __device__ void restoreU(
	const int &outerNormalX, const int &outerNormalY, const int &outerNormalZ,
	const float &rho, float &ux, float &uy, float &uz),
	const float (&f)[27]
)
{
	if ( outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[1] + f[3] + f[4] + f[5] + f[6] + (2.f) * f[7] + (2.f) * f[9] + (2.f) * f[12] + f[13] + f[14] + (2.f) * f[16] + f[17] + f[18] + (2.f) * f[20] + (2.f) * f[22] + (2.f) * f[23] + (2.f) * f[26];
		const float scmf1 = - f[3] + f[4] + (-2.f) * f[7] + (2.f) * f[9] - f[13] + f[14] + f[17] - f[18] + (2.f) * f[20] + (-2.f) * f[22] + (-2.f) * f[23] + (2.f) * f[26];
		const float scmf2 = - f[5] + f[6] + (2.f) * f[12] + f[13] - f[14] + (-2.f) * f[16] + f[17] - f[18] + (-2.f) * f[20] + (2.f) * f[22] + (-2.f) * f[23] + (2.f) * f[26];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (0.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0;
		uy = + s2;
		uz = + s1;
	}
	else if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + f[1] + f[2] + f[3] + f[4] + (2.f) * f[6] + f[7] + f[8] + f[9] + f[10] + (2.f) * f[12] + (2.f) * f[13] + (2.f) * f[15] + (2.f) * f[17] + (2.f) * f[19] + (2.f) * f[22] + (2.f) * f[24] + (2.f) * f[26];
		const float scmf1 = - f[3] + f[4] - f[7] + f[8] + f[9] - f[10] + (-2.f) * f[13] + (2.f) * f[17] + (-2.f) * f[19] + (-2.f) * f[22] + (2.f) * f[24] + (2.f) * f[26];
		const float scmf2 = + f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + (2.f) * f[12] + (-2.f) * f[15] + (-2.f) * f[19] + (2.f) * f[22] + (-2.f) * f[24] + (2.f) * f[26];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (0.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s2;
		uy = + s0;
		uz = + s1;
	}
	else if ( outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + f[1] + f[2] + (2.f) * f[4] + f[5] + f[6] + (2.f) * f[8] + (2.f) * f[9] + f[11] + f[12] + (2.f) * f[14] + f[15] + f[16] + (2.f) * f[17] + (2.f) * f[20] + (2.f) * f[21] + (2.f) * f[24] + (2.f) * f[26];
		const float scmf1 = - f[5] + f[6] - f[11] + f[12] + (-2.f) * f[14] + f[15] - f[16] + (2.f) * f[17] + (-2.f) * f[20] + (-2.f) * f[21] + (2.f) * f[24] + (2.f) * f[26];
		const float scmf2 = + f[1] - f[2] + (-2.f) * f[8] + (2.f) * f[9] - f[11] + f[12] - f[15] + f[16] + (2.f) * f[20] + (-2.f) * f[21] + (-2.f) * f[24] + (2.f) * f[26];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (0.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s2;
		uy = + s1;
		uz = + s0;
	}
	else if ( outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[2] - f[3] - f[4] - f[5] - f[6] + (-2.f) * f[8] + (-2.f) * f[10] + (-2.f) * f[11] - f[13] - f[14] + (-2.f) * f[15] - f[17] - f[18] + (-2.f) * f[19] + (-2.f) * f[21] + (-2.f) * f[24] + (-2.f) * f[25];
		const float scmf1 = + f[3] - f[4] + (-2.f) * f[8] + (2.f) * f[10] + f[13] - f[14] - f[17] + f[18] + (2.f) * f[19] + (-2.f) * f[21] + (-2.f) * f[24] + (2.f) * f[25];
		const float scmf2 = + f[5] - f[6] + (2.f) * f[11] - f[13] + f[14] + (-2.f) * f[15] - f[17] + f[18] + (-2.f) * f[19] + (2.f) * f[21] + (-2.f) * f[24] + (2.f) * f[25];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (0.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0;
		uy = - s2;
		uz = - s1;
	}
	else if ( outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] - f[1] - f[2] - f[3] - f[4] + (-2.f) * f[5] - f[7] - f[8] - f[9] - f[10] + (-2.f) * f[11] + (-2.f) * f[14] + (-2.f) * f[16] + (-2.f) * f[18] + (-2.f) * f[20] + (-2.f) * f[21] + (-2.f) * f[23] + (-2.f) * f[25];
		const float scmf1 = + f[3] - f[4] + f[7] - f[8] - f[9] + f[10] + (-2.f) * f[14] + (2.f) * f[18] + (-2.f) * f[20] + (-2.f) * f[21] + (2.f) * f[23] + (2.f) * f[25];
		const float scmf2 = - f[1] + f[2] - f[7] + f[8] - f[9] + f[10] + (2.f) * f[11] + (-2.f) * f[16] + (-2.f) * f[20] + (2.f) * f[21] + (-2.f) * f[23] + (2.f) * f[25];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (0.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s2;
		uy = + s0;
		uz = - s1;
	}
	else if ( outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] - f[1] - f[2] + (-2.f) * f[3] - f[5] - f[6] + (-2.f) * f[7] + (-2.f) * f[10] - f[11] - f[12] + (-2.f) * f[13] - f[15] - f[16] + (-2.f) * f[18] + (-2.f) * f[19] + (-2.f) * f[22] + (-2.f) * f[23] + (-2.f) * f[25];
		const float scmf1 = + f[5] - f[6] + f[11] - f[12] + (-2.f) * f[13] - f[15] + f[16] + (2.f) * f[18] + (-2.f) * f[19] + (-2.f) * f[22] + (2.f) * f[23] + (2.f) * f[25];
		const float scmf2 = - f[1] + f[2] + (-2.f) * f[7] + (2.f) * f[10] + f[11] - f[12] + f[15] - f[16] + (2.f) * f[19] + (-2.f) * f[22] + (-2.f) * f[23] + (2.f) * f[25];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (0.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s2;
		uy = - s1;
		uz = + s0;
	}
	else if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + f[1] + f[2] + (2.f) * f[4] + (2.f) * f[6] + (2.f) * f[8] + (2.f) * f[9] + (2.f) * f[12] + (2.f) * f[15] + (4.f) * f[17] + (4.f) * f[24] + (4.f) * f[26];
		const float scmf1 = - f[0] - f[1] - f[2] + (-2.f) * f[4] + (-2.f) * f[8] + (-2.f) * f[9];
		const float scmf2 = + f[1] - f[2] + (-2.f) * f[8] + (2.f) * f[9] + (2.f) * f[12] + (-2.f) * f[15] + (-4.f) * f[24] + (4.f) * f[26];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s2;
		uy = + s0 + (3.f/2.f) * s1;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[1] + (2.f) * f[4] + f[5] + f[6] + (4.f) * f[9] + (2.f) * f[12] + (2.f) * f[14] + (2.f) * f[16] + (2.f) * f[17] + (4.f) * f[20] + (4.f) * f[26];
		const float scmf1 = - f[0] + (-2.f) * f[4] - f[5] - f[6] + (-2.f) * f[14] + (-2.f) * f[17];
		const float scmf2 = - f[5] + f[6] + (2.f) * f[12] + (-2.f) * f[14] + (-2.f) * f[16] + (2.f) * f[17] + (-4.f) * f[20] + (4.f) * f[26];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (3.f/2.f) * s1;
		uy = + s2;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[1] + f[3] + f[4] + (2.f) * f[6] + (2.f) * f[7] + (2.f) * f[9] + (4.f) * f[12] + (2.f) * f[13] + (2.f) * f[17] + (4.f) * f[22] + (4.f) * f[26];
		const float scmf1 = - f[0] - f[3] - f[4] + (-2.f) * f[6] + (-2.f) * f[13] + (-2.f) * f[17];
		const float scmf2 = - f[3] + f[4] + (-2.f) * f[7] + (2.f) * f[9] + (-2.f) * f[13] + (2.f) * f[17] + (-4.f) * f[22] + (4.f) * f[26];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (3.f/2.f) * s1;
		uy = + (-3.f/2.f) * s1;
		uz = + s2;
	}
	else if ( outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + f[1] + f[2] + (2.f) * f[3] + (2.f) * f[5] + (2.f) * f[7] + (2.f) * f[10] + (2.f) * f[11] + (2.f) * f[16] + (4.f) * f[18] + (4.f) * f[23] + (4.f) * f[25];
		const float scmf1 = + f[0] + f[1] + f[2] + (2.f) * f[3] + (2.f) * f[7] + (2.f) * f[10];
		const float scmf2 = + f[1] - f[2] + (2.f) * f[7] + (-2.f) * f[10] + (-2.f) * f[11] + (2.f) * f[16] + (4.f) * f[23] + (-4.f) * f[25];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s2;
		uy = - s0 + (3.f/2.f) * s1;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[2] + (2.f) * f[3] + f[5] + f[6] + (4.f) * f[10] + (2.f) * f[11] + (2.f) * f[13] + (2.f) * f[15] + (2.f) * f[18] + (4.f) * f[19] + (4.f) * f[25];
		const float scmf1 = + f[0] + (2.f) * f[3] + f[5] + f[6] + (2.f) * f[13] + (2.f) * f[18];
		const float scmf2 = - f[5] + f[6] + (-2.f) * f[11] + (2.f) * f[13] + (2.f) * f[15] + (-2.f) * f[18] + (4.f) * f[19] + (-4.f) * f[25];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (3.f/2.f) * s1;
		uy = + s2;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[2] + f[3] + f[4] + (2.f) * f[5] + (2.f) * f[8] + (2.f) * f[10] + (4.f) * f[11] + (2.f) * f[14] + (2.f) * f[18] + (4.f) * f[21] + (4.f) * f[25];
		const float scmf1 = + f[0] + f[3] + f[4] + (2.f) * f[5] + (2.f) * f[14] + (2.f) * f[18];
		const float scmf2 = - f[3] + f[4] + (2.f) * f[8] + (-2.f) * f[10] + (2.f) * f[14] + (-2.f) * f[18] + (4.f) * f[21] + (-4.f) * f[25];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (3.f/2.f) * s1;
		uy = + (-3.f/2.f) * s1;
		uz = + s2;
	}
	else if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] - f[1] - f[2] + (-2.f) * f[3] + (-2.f) * f[6] + (-2.f) * f[7] + (-2.f) * f[10] + (-2.f) * f[12] + (-4.f) * f[13] + (-2.f) * f[15] + (-4.f) * f[19] + (-4.f) * f[22];
		const float scmf1 = + f[0] + f[1] + f[2] + (2.f) * f[3] + (2.f) * f[7] + (2.f) * f[10];
		const float scmf2 = - f[1] + f[2] + (-2.f) * f[7] + (2.f) * f[10] + (-2.f) * f[12] + (2.f) * f[15] + (4.f) * f[19] + (-4.f) * f[22];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s2;
		uy = - s0 + (-3.f/2.f) * s1;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[1] + (-2.f) * f[3] - f[5] - f[6] + (-4.f) * f[7] + (-2.f) * f[12] + (-2.f) * f[13] + (-2.f) * f[16] + (-2.f) * f[18] + (-4.f) * f[22] + (-4.f) * f[23];
		const float scmf1 = + f[0] + (2.f) * f[3] + f[5] + f[6] + (2.f) * f[13] + (2.f) * f[18];
		const float scmf2 = + f[5] - f[6] + (-2.f) * f[12] + (-2.f) * f[13] + (2.f) * f[16] + (2.f) * f[18] + (-4.f) * f[22] + (4.f) * f[23];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (-3.f/2.f) * s1;
		uy = - s2;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[1] - f[3] - f[4] + (-2.f) * f[5] + (-2.f) * f[7] + (-2.f) * f[9] + (-2.f) * f[14] + (-4.f) * f[16] + (-2.f) * f[18] + (-4.f) * f[20] + (-4.f) * f[23];
		const float scmf1 = + f[0] + f[3] + f[4] + (2.f) * f[5] + (2.f) * f[14] + (2.f) * f[18];
		const float scmf2 = + f[3] - f[4] + (2.f) * f[7] + (-2.f) * f[9] + (-2.f) * f[14] + (2.f) * f[18] + (-4.f) * f[20] + (4.f) * f[23];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (-3.f/2.f) * s1;
		uy = + (-3.f/2.f) * s1;
		uz = - s2;
	}
	else if ( outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] - f[1] - f[2] + (-2.f) * f[4] + (-2.f) * f[5] + (-2.f) * f[8] + (-2.f) * f[9] + (-2.f) * f[11] + (-4.f) * f[14] + (-2.f) * f[16] + (-4.f) * f[20] + (-4.f) * f[21];
		const float scmf1 = - f[0] - f[1] - f[2] + (-2.f) * f[4] + (-2.f) * f[8] + (-2.f) * f[9];
		const float scmf2 = - f[1] + f[2] + (2.f) * f[8] + (-2.f) * f[9] + (2.f) * f[11] + (-2.f) * f[16] + (-4.f) * f[20] + (4.f) * f[21];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s2;
		uy = + s0 + (-3.f/2.f) * s1;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[2] + (-2.f) * f[4] - f[5] - f[6] + (-4.f) * f[8] + (-2.f) * f[11] + (-2.f) * f[14] + (-2.f) * f[15] + (-2.f) * f[17] + (-4.f) * f[21] + (-4.f) * f[24];
		const float scmf1 = - f[0] + (-2.f) * f[4] - f[5] - f[6] + (-2.f) * f[14] + (-2.f) * f[17];
		const float scmf2 = + f[5] - f[6] + (2.f) * f[11] + (2.f) * f[14] + (-2.f) * f[15] + (-2.f) * f[17] + (4.f) * f[21] + (-4.f) * f[24];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (-3.f/2.f) * s1;
		uy = - s2;
		uz = + (-3.f/2.f) * s1;
	}
	else if ( outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 0 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[2] - f[3] - f[4] + (-2.f) * f[6] + (-2.f) * f[8] + (-2.f) * f[10] + (-2.f) * f[13] + (-4.f) * f[15] + (-2.f) * f[17] + (-4.f) * f[19] + (-4.f) * f[24];
		const float scmf1 = - f[0] - f[3] - f[4] + (-2.f) * f[6] + (-2.f) * f[13] + (-2.f) * f[17];
		const float scmf2 = + f[3] - f[4] + (-2.f) * f[8] + (2.f) * f[10] + (2.f) * f[13] + (-2.f) * f[17] + (4.f) * f[19] + (-4.f) * f[24];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (0.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (-3.f/2.f) * s1;
		uy = + (-3.f/2.f) * s1;
		uz = - s2;
	}
	else if ( outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[1] + (2.f) * f[4] + (2.f) * f[6] + (4.f) * f[9] + (4.f) * f[12] + (4.f) * f[17] + (8.f) * f[26];
		const float scmf1 = - f[0] + (-2.f) * f[4] + (-2.f) * f[6] + (-4.f) * f[17];
		const float scmf2 = - f[0] + (-2.f) * f[1] + (-2.f) * f[4] + (-4.f) * f[9];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (-2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (3.f/2.f) * s1;
		uy = + s0 + (3.f/2.f) * s2;
		uz = - s0 + (-3.f/2.f) * s1 + (-3.f/2.f) * s2;
	}
	else if ( outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[2] + (-2.f) * f[3] + (-2.f) * f[5] + (-4.f) * f[10] + (-4.f) * f[11] + (-4.f) * f[18] + (-8.f) * f[25];
		const float scmf1 = - f[0] + (-2.f) * f[3] + (-2.f) * f[5] + (-4.f) * f[18];
		const float scmf2 = - f[0] + (-2.f) * f[2] + (-2.f) * f[3] + (-4.f) * f[10];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (-2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (-3.f/2.f) * s1;
		uy = + s0 + (-3.f/2.f) * s2;
		uz = - s0 + (3.f/2.f) * s1 + (3.f/2.f) * s2;
	}
	else if ( outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[2] + (-2.f) * f[4] + (-2.f) * f[6] + (-4.f) * f[8] + (-4.f) * f[15] + (-4.f) * f[17] + (-8.f) * f[24];
		const float scmf1 = - f[0] + (-2.f) * f[4] + (-2.f) * f[6] + (-4.f) * f[17];
		const float scmf2 = + f[0] + (2.f) * f[2] + (2.f) * f[4] + (4.f) * f[8];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (-3.f/2.f) * s1;
		uy = - s0 + (-3.f/2.f) * s2;
		uz = + s0 + (-3.f/2.f) * s1 + (3.f/2.f) * s2;
	}
	else if ( outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[1] + (-2.f) * f[4] + (-2.f) * f[5] + (-4.f) * f[9] + (-4.f) * f[14] + (-4.f) * f[16] + (-8.f) * f[20];
		const float scmf1 = + f[0] + (2.f) * f[4] + (2.f) * f[5] + (4.f) * f[14];
		const float scmf2 = - f[0] + (-2.f) * f[1] + (-2.f) * f[4] + (-4.f) * f[9];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (-2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (-3.f/2.f) * s1;
		uy = + s0 + (-3.f/2.f) * s2;
		uz = + s0 + (3.f/2.f) * s1 + (-3.f/2.f) * s2;
	}
	else if ( outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = - f[0] + (-2.f) * f[1] + (-2.f) * f[3] + (-2.f) * f[6] + (-4.f) * f[7] + (-4.f) * f[12] + (-4.f) * f[13] + (-8.f) * f[22];
		const float scmf1 = + f[0] + (2.f) * f[3] + (2.f) * f[6] + (4.f) * f[13];
		const float scmf2 = + f[0] + (2.f) * f[1] + (2.f) * f[3] + (4.f) * f[7];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (-1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (-3.f/2.f) * s1;
		uy = - s0 + (-3.f/2.f) * s2;
		uz = - s0 + (-3.f/2.f) * s1 + (-3.f/2.f) * s2;
	}
	else if ( outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[1] + (2.f) * f[3] + (2.f) * f[5] + (4.f) * f[7] + (4.f) * f[16] + (4.f) * f[18] + (8.f) * f[23];
		const float scmf1 = - f[0] + (-2.f) * f[3] + (-2.f) * f[5] + (-4.f) * f[18];
		const float scmf2 = + f[0] + (2.f) * f[1] + (2.f) * f[3] + (4.f) * f[7];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (-2.f/3.f);
		const float s2 = scmf2 - (2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = + s0 + (3.f/2.f) * s1;
		uy = - s0 + (3.f/2.f) * s2;
		uz = + s0 + (3.f/2.f) * s1 + (-3.f/2.f) * s2;
	}
	else if ( outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == -1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[2] + (2.f) * f[3] + (2.f) * f[6] + (4.f) * f[10] + (4.f) * f[13] + (4.f) * f[15] + (8.f) * f[19];
		const float scmf1 = + f[0] + (2.f) * f[3] + (2.f) * f[6] + (4.f) * f[13];
		const float scmf2 = - f[0] + (-2.f) * f[2] + (-2.f) * f[3] + (-4.f) * f[10];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (-2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (3.f/2.f) * s1;
		uy = + s0 + (3.f/2.f) * s2;
		uz = + s0 + (-3.f/2.f) * s1 + (3.f/2.f) * s2;
	}
	else if ( outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 1 )
	{
		// Multiply SCMfk fk
		const float scmf0 = + f[0] + (2.f) * f[2] + (2.f) * f[4] + (2.f) * f[5] + (4.f) * f[8] + (4.f) * f[11] + (4.f) * f[14] + (8.f) * f[21];
		const float scmf1 = + f[0] + (2.f) * f[4] + (2.f) * f[5] + (4.f) * f[14];
		const float scmf2 = + f[0] + (2.f) * f[2] + (2.f) * f[4] + (4.f) * f[8];
		// Divide by rho
		scmf0 = scmf0 / rho;
		scmf1 = scmf1 / rho;
		scmf2 = scmf2 / rho;
		// Subtract 1/rho scmf - scq
		const float s0 = scmf0 - (1.f);
		const float s1 = scmf1 - (2.f/3.f);
		const float s2 = scmf2 - (2.f/3.f);
		// Multiply (Su C^T Qu)^-1 * s to get u
		ux = - s0 + (3.f/2.f) * s1;
		uy = - s0 + (3.f/2.f) * s2;
		uz = - s0 + (3.f/2.f) * s1 + (3.f/2.f) * s2;
	}
}