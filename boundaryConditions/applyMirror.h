__host__ __device__ void applyMirror(
	const int &outerNormalX, const int &outerNormalY, const int &outerNormalZ,
	float (&f)[27]
)
{
	if ( outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 0 )
	{	
	f[2] = f[1];
	f[8] = f[9];
	f[10] = f[7];
	f[11] = f[16];
	f[15] = f[12];
	f[19] = f[22];
	f[21] = f[20];
	f[24] = f[26];
	f[25] = f[23];
	}
	if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0 )
	{
	f[5] = f[6];
	f[11] = f[15];
	f[14] = f[17];
	f[16] = f[12];
	f[18] = f[13];
	f[20] = f[26];
	f[21] = f[24];
	f[23] = f[22];
	f[25] = f[19];
	}
}
