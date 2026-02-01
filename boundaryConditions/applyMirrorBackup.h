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
	else if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 0 )
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
	else if ( outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == 1 )
	{
	f[3] = f[4];
	f[7] = f[9];
	f[10] = f[8];
	f[13] = f[17];
	f[18] = f[14];
	f[19] = f[24];
	f[22] = f[26];
	f[23] = f[20];
	f[25] = f[21];
	}
	else if ( outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 0 )
	{
	f[1] = f[2];
	f[7] = f[10];
	f[9] = f[8];
	f[12] = f[15];
	f[16] = f[11];
	f[20] = f[21];
	f[22] = f[19];
	f[23] = f[25];
	f[26] = f[24];
	}
	else if ( outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 0 )
	{
	f[6] = f[5];
	f[12] = f[16];
	f[13] = f[18];
	f[15] = f[11];
	f[17] = f[14];
	f[19] = f[25];
	f[22] = f[23];
	f[24] = f[21];
	f[26] = f[20];
	}
	else if ( outerNormalX == 0 && outerNormalY == 0 && outerNormalZ == -1 )
	{
	f[4] = f[3];
	f[8] = f[10];
	f[9] = f[7];
	f[14] = f[18];
	f[17] = f[13];
	f[20] = f[23];
	f[21] = f[25];
	f[24] = f[19];
	f[26] = f[22];
	}
	else if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == 1 )
	{
	f[3] = f[4];
	f[5] = f[6];
	f[7] = f[9];
	f[10] = f[8];
	f[11] = f[15];
	f[13] = f[17];
	f[14] = f[17];
	f[16] = f[12];
	f[18] = f[17];
	f[19] = f[24];
	f[20] = f[26];
	f[21] = f[24];
	f[22] = f[26];
	f[23] = f[26];
	f[25] = f[24];
	}
	else if ( outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == 1 )
	{
	f[2] = f[1];
	f[3] = f[4];
	f[7] = f[9];
	f[8] = f[9];
	f[10] = f[9];
	f[11] = f[16];
	f[13] = f[17];
	f[15] = f[12];
	f[18] = f[14];
	f[19] = f[26];
	f[21] = f[20];
	f[22] = f[26];
	f[23] = f[20];
	f[24] = f[26];
	f[25] = f[20];
	}
	else if ( outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 0 )
	{
	f[2] = f[1];
	f[5] = f[6];
	f[8] = f[9];
	f[10] = f[7];
	f[11] = f[12];
	f[14] = f[17];
	f[15] = f[12];
	f[16] = f[12];
	f[18] = f[13];
	f[19] = f[22];
	f[20] = f[26];
	f[21] = f[26];
	f[23] = f[22];
	f[24] = f[26];
	f[25] = f[22];
	}
	else if ( outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == -1 )
	{
	f[4] = f[3];
	f[6] = f[5];
	f[8] = f[10];
	f[9] = f[7];
	f[12] = f[16];
	f[13] = f[18];
	f[14] = f[18];
	f[15] = f[11];
	f[17] = f[18];
	f[19] = f[25];
	f[20] = f[23];
	f[21] = f[25];
	f[22] = f[23];
	f[24] = f[25];
	f[26] = f[23];
	}
	else if ( outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == -1 )
	{
	f[1] = f[2];
	f[4] = f[3];
	f[7] = f[10];
	f[8] = f[10];
	f[9] = f[10];
	f[12] = f[15];
	f[14] = f[18];
	f[16] = f[11];
	f[17] = f[13];
	f[20] = f[25];
	f[21] = f[25];
	f[22] = f[19];
	f[23] = f[25];
	f[24] = f[19];
	f[26] = f[19];
	}
	else if ( outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 0 )
	{
	f[1] = f[2];
	f[6] = f[5];
	f[7] = f[10];
	f[9] = f[8];
	f[12] = f[11];
	f[13] = f[18];
	f[15] = f[11];
	f[16] = f[11];
	f[17] = f[14];
	f[19] = f[25];
	f[20] = f[21];
	f[22] = f[25];
	f[23] = f[25];
	f[24] = f[21];
	f[26] = f[21];
	}
	else if ( outerNormalX == 0 && outerNormalY == 1 && outerNormalZ == -1 )
	{
	f[4] = f[3];
	f[5] = f[6];
	f[8] = f[10];
	f[9] = f[7];
	f[11] = f[15];
	f[14] = f[13];
	f[16] = f[12];
	f[17] = f[13];
	f[18] = f[13];
	f[20] = f[22];
	f[21] = f[19];
	f[23] = f[22];
	f[24] = f[19];
	f[25] = f[19];
	f[26] = f[22];
	}
	else if ( outerNormalX == 1 && outerNormalY == 0 && outerNormalZ == -1 )
	{
	f[2] = f[1];
	f[4] = f[3];
	f[8] = f[7];
	f[9] = f[7];
	f[10] = f[7];
	f[11] = f[16];
	f[14] = f[18];
	f[15] = f[12];
	f[17] = f[13];
	f[19] = f[22];
	f[20] = f[23];
	f[21] = f[23];
	f[24] = f[22];
	f[25] = f[23];
	f[26] = f[22];
	}
	else if ( outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 0 )
	{
	f[2] = f[1];
	f[6] = f[5];
	f[8] = f[9];
	f[10] = f[7];
	f[11] = f[16];
	f[12] = f[16];
	f[13] = f[18];
	f[15] = f[16];
	f[17] = f[14];
	f[19] = f[23];
	f[21] = f[20];
	f[22] = f[23];
	f[24] = f[20];
	f[25] = f[23];
	f[26] = f[20];
	}
	else if ( outerNormalX == 0 && outerNormalY == -1 && outerNormalZ == 1 )
	{
	f[3] = f[4];
	f[6] = f[5];
	f[7] = f[9];
	f[10] = f[8];
	f[12] = f[16];
	f[13] = f[14];
	f[15] = f[11];
	f[17] = f[14];
	f[18] = f[14];
	f[19] = f[21];
	f[22] = f[20];
	f[23] = f[20];
	f[24] = f[21];
	f[25] = f[21];
	f[26] = f[20];
	}
	else if ( outerNormalX == -1 && outerNormalY == 0 && outerNormalZ == 1 )
	{
	f[1] = f[2];
	f[3] = f[4];
	f[7] = f[8];
	f[9] = f[8];
	f[10] = f[8];
	f[12] = f[15];
	f[13] = f[17];
	f[16] = f[11];
	f[18] = f[14];
	f[19] = f[24];
	f[20] = f[21];
	f[22] = f[24];
	f[23] = f[21];
	f[25] = f[21];
	f[26] = f[24];
	}
	else if ( outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 0 )
	{
	f[1] = f[2];
	f[5] = f[6];
	f[7] = f[10];
	f[9] = f[8];
	f[11] = f[15];
	f[12] = f[15];
	f[14] = f[17];
	f[16] = f[15];
	f[18] = f[13];
	f[20] = f[24];
	f[21] = f[24];
	f[22] = f[19];
	f[23] = f[19];
	f[25] = f[19];
	f[26] = f[24];
	}
	else if ( outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == 1 )
	{
	f[2] = f[1];
	f[3] = f[4];
	f[5] = f[6];
	f[7] = f[9];
	f[8] = f[9];
	f[10] = f[9];
	f[11] = f[12];
	f[13] = f[17];
	f[14] = f[17];
	f[15] = f[12];
	f[16] = f[12];
	f[18] = f[17];
	f[19] = f[26];
	f[20] = f[26];
	f[21] = f[26];
	f[22] = f[26];
	f[23] = f[26];
	f[24] = f[26];
	f[25] = f[26];
	}
	else if ( outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == -1 )
	{
	f[1] = f[2];
	f[4] = f[3];
	f[6] = f[5];
	f[7] = f[10];
	f[8] = f[10];
	f[9] = f[10];
	f[12] = f[11];
	f[13] = f[18];
	f[14] = f[18];
	f[15] = f[11];
	f[16] = f[11];
	f[17] = f[18];
	f[19] = f[25];
	f[20] = f[25];
	f[21] = f[25];
	f[22] = f[25];
	f[23] = f[25];
	f[24] = f[25];
	f[26] = f[25];
	}
	else if ( outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == 1 )
	{
	f[1] = f[2];
	f[3] = f[4];
	f[5] = f[6];
	f[7] = f[8];
	f[9] = f[8];
	f[10] = f[8];
	f[11] = f[15];
	f[12] = f[15];
	f[13] = f[17];
	f[14] = f[17];
	f[16] = f[15];
	f[18] = f[17];
	f[19] = f[24];
	f[20] = f[24];
	f[21] = f[24];
	f[22] = f[24];
	f[23] = f[24];
	f[25] = f[24];
	f[26] = f[24];
	}
	else if ( outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == 1 )
	{
	f[2] = f[1];
	f[3] = f[4];
	f[6] = f[5];
	f[7] = f[9];
	f[8] = f[9];
	f[10] = f[9];
	f[11] = f[16];
	f[12] = f[16];
	f[13] = f[14];
	f[15] = f[16];
	f[17] = f[14];
	f[18] = f[14];
	f[19] = f[20];
	f[21] = f[20];
	f[22] = f[20];
	f[23] = f[20];
	f[24] = f[20];
	f[25] = f[20];
	f[26] = f[20];
	}
	else if ( outerNormalX == 1 && outerNormalY == 1 && outerNormalZ == -1 )
	{
	f[2] = f[1];
	f[4] = f[3];
	f[5] = f[6];
	f[8] = f[7];
	f[9] = f[7];
	f[10] = f[7];
	f[11] = f[12];
	f[14] = f[13];
	f[15] = f[12];
	f[16] = f[12];
	f[17] = f[13];
	f[18] = f[13];
	f[19] = f[22];
	f[20] = f[22];
	f[21] = f[22];
	f[23] = f[22];
	f[24] = f[22];
	f[25] = f[22];
	f[26] = f[22];
	}
	else if ( outerNormalX == 1 && outerNormalY == -1 && outerNormalZ == -1 )
	{
	f[2] = f[1];
	f[4] = f[3];
	f[6] = f[5];
	f[8] = f[7];
	f[9] = f[7];
	f[10] = f[7];
	f[11] = f[16];
	f[12] = f[16];
	f[13] = f[18];
	f[14] = f[18];
	f[15] = f[16];
	f[17] = f[18];
	f[19] = f[23];
	f[20] = f[23];
	f[21] = f[23];
	f[22] = f[23];
	f[24] = f[23];
	f[25] = f[23];
	f[26] = f[23];
	}
	else if ( outerNormalX == -1 && outerNormalY == 1 && outerNormalZ == -1 )
	{
	f[1] = f[2];
	f[4] = f[3];
	f[5] = f[6];
	f[7] = f[10];
	f[8] = f[10];
	f[9] = f[10];
	f[11] = f[15];
	f[12] = f[15];
	f[14] = f[13];
	f[16] = f[15];
	f[17] = f[13];
	f[18] = f[13];
	f[20] = f[19];
	f[21] = f[19];
	f[22] = f[19];
	f[23] = f[19];
	f[24] = f[19];
	f[25] = f[19];
	f[26] = f[19];
	}
	else if ( outerNormalX == -1 && outerNormalY == -1 && outerNormalZ == 1 )
	{
	f[1] = f[2];
	f[3] = f[4];
	f[6] = f[5];
	f[7] = f[8];
	f[9] = f[8];
	f[10] = f[8];
	f[12] = f[11];
	f[13] = f[14];
	f[15] = f[11];
	f[16] = f[11];
	f[17] = f[14];
	f[18] = f[14];
	f[19] = f[21];
	f[20] = f[21];
	f[22] = f[21];
	f[23] = f[21];
	f[24] = f[21];
	f[25] = f[21];
	f[26] = f[21];
	}
}
