#pragma once

// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

/// Raoyang Zhang, Hongli Fan, Hudong Chen - A lattice Boltzmann approach for solving scalar transport equations, 2011
__host__ __device__ void applyScalarTransportCollision( float (&f)[27], float (&T)[27], const float &tauT )
{
	float rho, ux, uy, uz;
	getRhoUxUyUz( rho, ux, uy, uz, f );
	
	// eq 2.7
	float Tsum = 0.f;
	for ( int direction = 0; direction < 27; direction++ ) Tsum += f[direction] * T[direction];
	Tsum = Tsum / rho;
	
	// eq 2.6
	float fjTjminT[27];
	for ( int direction = 0; direction < 27; direction++ ) fjTjminT[direction] = f[direction] * (T[direction] - Tsum);
	
	const float sumX = 	+fjTjminT[1]  - fjTjminT[2]  + fjTjminT[7]  - fjTjminT[8]  + fjTjminT[9]  - fjTjminT[10]
						-fjTjminT[11] + fjTjminT[12] - fjTjminT[15] + fjTjminT[16] - fjTjminT[19] + fjTjminT[20]
						-fjTjminT[21] + fjTjminT[22] + fjTjminT[23] - fjTjminT[24] - fjTjminT[25] + fjTjminT[26];

	const float sumY = 	-fjTjminT[5]  + fjTjminT[6]  - fjTjminT[11] + fjTjminT[12] + fjTjminT[13] - fjTjminT[14]
						+fjTjminT[15] - fjTjminT[16] + fjTjminT[17] - fjTjminT[18] + fjTjminT[19] - fjTjminT[20]
						-fjTjminT[21] + fjTjminT[22] - fjTjminT[23] + fjTjminT[24] - fjTjminT[25] + fjTjminT[26];

	const float sumZ = 	-fjTjminT[3]  + fjTjminT[4]  - fjTjminT[7]  + fjTjminT[8]  + fjTjminT[9]  - fjTjminT[10]
						-fjTjminT[13] + fjTjminT[14] + fjTjminT[17] - fjTjminT[18] - fjTjminT[19] + fjTjminT[20]
						+fjTjminT[21] - fjTjminT[22] - fjTjminT[23] + fjTjminT[24] - fjTjminT[25] + fjTjminT[26];
	
	const float rhodiv = 3.f / rho;
	
	float Fi[27];
	Fi[0]  = rhodiv * ( (     - ux)*sumX + (     - uy)*sumY + (     - uz)*sumZ );
	Fi[1]  = rhodiv * ( ( 1.f - ux)*sumX + (     - uy)*sumY + (     - uz)*sumZ );
	Fi[2]  = rhodiv * ( (-1.f - ux)*sumX + (     - uy)*sumY + (     - uz)*sumZ );
	Fi[3]  = rhodiv * ( (     - ux)*sumX + (     - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[4]  = rhodiv * ( (     - ux)*sumX + (     - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[5]  = rhodiv * ( (     - ux)*sumX + (-1.f - uy)*sumY + (     - uz)*sumZ );
	Fi[6]  = rhodiv * ( (     - ux)*sumX + ( 1.f - uy)*sumY + (     - uz)*sumZ );
	Fi[7]  = rhodiv * ( ( 1.f - ux)*sumX + (     - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[8]  = rhodiv * ( (-1.f - ux)*sumX + (     - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[9]  = rhodiv * ( ( 1.f - ux)*sumX + (     - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[10] = rhodiv * ( (-1.f - ux)*sumX + (     - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[11] = rhodiv * ( (-1.f - ux)*sumX + (-1.f - uy)*sumY + (     - uz)*sumZ );
	Fi[12] = rhodiv * ( ( 1.f - ux)*sumX + ( 1.f - uy)*sumY + (     - uz)*sumZ );
	Fi[13] = rhodiv * ( (     - ux)*sumX + ( 1.f - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[14] = rhodiv * ( (     - ux)*sumX + (-1.f - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[15] = rhodiv * ( (-1.f - ux)*sumX + ( 1.f - uy)*sumY + (     - uz)*sumZ );
	Fi[16] = rhodiv * ( ( 1.f - ux)*sumX + (-1.f - uy)*sumY + (     - uz)*sumZ );
	Fi[17] = rhodiv * ( (     - ux)*sumX + ( 1.f - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[18] = rhodiv * ( (     - ux)*sumX + (-1.f - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[19] = rhodiv * ( (-1.f - ux)*sumX + ( 1.f - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[20] = rhodiv * ( ( 1.f - ux)*sumX + (-1.f - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[21] = rhodiv * ( (-1.f - ux)*sumX + (-1.f - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[22] = rhodiv * ( ( 1.f - ux)*sumX + ( 1.f - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[23] = rhodiv * ( ( 1.f - ux)*sumX + (-1.f - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[24] = rhodiv * ( (-1.f - ux)*sumX + ( 1.f - uy)*sumY + ( 1.f - uz)*sumZ );
	Fi[25] = rhodiv * ( (-1.f - ux)*sumX + (-1.f - uy)*sumY + (-1.f - uz)*sumZ );
	Fi[26] = rhodiv * ( ( 1.f - ux)*sumX + ( 1.f - uy)*sumY + ( 1.f - uz)*sumZ );
	
	// eq 2.5
	const float bracket = 1.f - 1.f / tauT;
	for ( int direction = 0; direction < 27; direction++ ) T[direction] = Tsum + bracket * Fi[direction];
}

__host__ __device__ void applyScalarTransportGivenT( float (&T)[27], float &givenT )
{
	for ( int direction = 0; direction < 27; direction++ ) T[direction] = givenT;
}
