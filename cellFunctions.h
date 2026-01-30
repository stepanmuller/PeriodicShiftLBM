#pragma once

// id: { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: { 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: { 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

// cx * cx: { 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 };
// cy * cy: { 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
// cz * cz: { 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

// cy * cz: { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1 };
// cx * cz: { 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1,-1,-1, 1, 1 };
// cx * cy: { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,-1,-1, 0, 0,-1,-1, 1, 1,-1,-1, 1, 1 };

// w:  { 8/27, 2/27, 2/27, 2/27 , 2/27, 2/27, 2/27, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/54, 1/216, 1/216, 1/216, 1/216, 1/216, 1/216, 1/216, 1/216 };

__host__ __device__ void getFeq(
	const float &rho, const float &ux, const float &uy, const float &uz, 
	float (&feq)[27]
	)
{

	const float u2 = ux*ux + uy*uy + uz*uz;

	const float cu0  = 0.f;
	const float cu1  = +ux;
	const float cu2  = -ux;
	const float cu3  = -uz;
	const float cu4  = +uz;
	const float cu5  = -uy;
	const float cu6  = +uy;
	const float cu7  = +ux -uz;
	const float cu8  = -ux +uz;
	const float cu9  = +ux +uz;
	const float cu10 = -ux -uz;
	const float cu11 = -ux -uy;
	const float cu12 = +ux +uy;
	const float cu13 = +uy -uz;
	const float cu14 = -uy +uz;
	const float cu15 = -ux +uy;
	const float cu16 = +ux -uy;
	const float cu17 = +uy +uz;
	const float cu18 = -uy -uz;
	const float cu19 = -ux +uy -uz;
	const float cu20 = +ux -uy +uz;
	const float cu21 = -ux -uy +uz;
	const float cu22 = +ux +uy -uz;
	const float cu23 = +ux -uy -uz;
	const float cu24 = -ux +uy +uz;
	const float cu25 = -ux -uy -uz;
	const float cu26 = +ux +uy +uz;

	constexpr float w0  = 8.f/27.f;
	constexpr float w1  = 2.f/27.f;
	constexpr float w2  = 1.f/54.f;
	constexpr float w3 = 1.f/216.f;

	feq[0]  = rho*w0 *(1.f + 3.f*cu0  + 4.5f*cu0 *cu0  - 1.5f*u2);
	feq[1]  = rho*w1 *(1.f + 3.f*cu1  + 4.5f*cu1 *cu1  - 1.5f*u2);
	feq[2]  = rho*w1 *(1.f + 3.f*cu2  + 4.5f*cu2 *cu2  - 1.5f*u2);
	feq[3]  = rho*w1 *(1.f + 3.f*cu3  + 4.5f*cu3 *cu3  - 1.5f*u2);
	feq[4]  = rho*w1 *(1.f + 3.f*cu4  + 4.5f*cu4 *cu4  - 1.5f*u2);
	feq[5]  = rho*w1 *(1.f + 3.f*cu5  + 4.5f*cu5 *cu5  - 1.5f*u2);
	feq[6]  = rho*w1 *(1.f + 3.f*cu6  + 4.5f*cu6 *cu6  - 1.5f*u2);
	feq[7]  = rho*w2 *(1.f + 3.f*cu7  + 4.5f*cu7 *cu7  - 1.5f*u2);
	feq[8]  = rho*w2 *(1.f + 3.f*cu8  + 4.5f*cu8 *cu8  - 1.5f*u2);
	feq[9]  = rho*w2 *(1.f + 3.f*cu9  + 4.5f*cu9 *cu9  - 1.5f*u2);
	feq[10] = rho*w2 *(1.f + 3.f*cu10 + 4.5f*cu10*cu10 - 1.5f*u2);
	feq[11] = rho*w2 *(1.f + 3.f*cu11 + 4.5f*cu11*cu11 - 1.5f*u2);
	feq[12] = rho*w2 *(1.f + 3.f*cu12 + 4.5f*cu12*cu12 - 1.5f*u2);
	feq[13] = rho*w2 *(1.f + 3.f*cu13 + 4.5f*cu13*cu13 - 1.5f*u2);
	feq[14] = rho*w2 *(1.f + 3.f*cu14 + 4.5f*cu14*cu14 - 1.5f*u2);
	feq[15] = rho*w2 *(1.f + 3.f*cu15 + 4.5f*cu15*cu15 - 1.5f*u2);
	feq[16] = rho*w2 *(1.f + 3.f*cu16 + 4.5f*cu16*cu16 - 1.5f*u2);
	feq[17] = rho*w2 *(1.f + 3.f*cu17 + 4.5f*cu17*cu17 - 1.5f*u2);
	feq[18] = rho*w2 *(1.f + 3.f*cu18 + 4.5f*cu18*cu18 - 1.5f*u2);
	feq[19] = rho*w3 *(1.f + 3.f*cu19 + 4.5f*cu19*cu19 - 1.5f*u2);
	feq[20] = rho*w3 *(1.f + 3.f*cu20 + 4.5f*cu20*cu20 - 1.5f*u2);
	feq[21] = rho*w3 *(1.f + 3.f*cu21 + 4.5f*cu21*cu21 - 1.5f*u2);
	feq[22] = rho*w3 *(1.f + 3.f*cu22 + 4.5f*cu22*cu22 - 1.5f*u2);
	feq[23] = rho*w3 *(1.f + 3.f*cu23 + 4.5f*cu23*cu23 - 1.5f*u2);
	feq[24] = rho*w3 *(1.f + 3.f*cu24 + 4.5f*cu24*cu24 - 1.5f*u2);
	feq[25] = rho*w3 *(1.f + 3.f*cu25 + 4.5f*cu25*cu25 - 1.5f*u2);
	feq[26] = rho*w3 *(1.f + 3.f*cu26 + 4.5f*cu26*cu26 - 1.5f*u2);	
}

__host__ __device__ void getFneq(const float (&f)[27], const float (&feq)[27], float (&fneq)[27])
{
	for (size_t i = 0; i < 27; i++)	fneq[i] = f[i] - feq[i];
}

__host__ __device__ void getRhoUxUyUz(
	float &rho, float &ux, float &uy, float &uz, 
	const float (&f)[27]
	)
{
	rho = 	f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] + f[9] +
			f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18] + f[19] +
			f[20] + f[21] + f[22] + f[23] + f[24] + f[25] + f[26];     
	
	ux = 	(+f[1]  - f[2]  + f[7]  - f[8]  + f[9]  - f[10]
			-f[11] + f[12] - f[15] + f[16] - f[19] + f[20]
			-f[21] + f[22] + f[23] - f[24] - f[25] + f[26]) / rho;

	uy = 	(-f[5]  + f[6]  - f[11] + f[12] + f[13] - f[14]
			+f[15] - f[16] + f[17] - f[18] + f[19] - f[20]
			-f[21] + f[22] - f[23] + f[24] - f[25] + f[26]) / rho;

	uz = 	(-f[3]  + f[4]  - f[7]  + f[8]  + f[9]  - f[10]
			-f[13] + f[14] + f[17] - f[18] - f[19] + f[20]
			+f[21] - f[22] - f[23] + f[24] - f[25] + f[26]) / rho;
}

__host__ __device__ void getOmegaLES(const float (&fneq)[27], const float &rho, float &omegaLES)
{
	if (SmagorinskyConstant == 0)
	{
		omegaLES = 1 / tau;
		return;
	}
	
	float P = 0.f;

	const float cxcx = (fneq[1] + fneq[2] + fneq[7] + fneq[8] + fneq[9] + fneq[10] + fneq[11] + fneq[12] + fneq[15] 
				+ fneq[16] + fneq[19] + fneq[20] + fneq[21] + fneq[22] + fneq[23] + fneq[24] + fneq[25] + fneq[26]);
	P += (cxcx * cxcx);

	const float cycy = (fneq[5] + fneq[6] + fneq[11] + fneq[12] + fneq[13] + fneq[14] + fneq[15] + fneq[16] + fneq[17]
				+ fneq[18] + fneq[19] + fneq[20] + fneq[21] + fneq[22] + fneq[23] + fneq[24] + fneq[25] + fneq[26]);
	P += (cycy * cycy);

	const float czcz = (fneq[3] + fneq[4] + fneq[7] + fneq[8] + fneq[9] + fneq[10] + fneq[13] + fneq[14] + fneq[17]
				+ fneq[18] + fneq[19] + fneq[20] + fneq[21] + fneq[22] + fneq[23] + fneq[24] + fneq[25] + fneq[26]);
	P += (czcz * czcz);

	const float cycz = (-fneq[13] - fneq[14] + fneq[17] + fneq[18] - fneq[19] - fneq[20] - fneq[21] - fneq[22]
				+ fneq[23] + fneq[24] + fneq[25] + fneq[26]);
	P += 2.f*(cycz * cycz);

	const float cxcz = (-fneq[7] - fneq[8] + fneq[9] + fneq[10] + fneq[19] + fneq[20] - fneq[21] - fneq[22]
				- fneq[23] - fneq[24] + fneq[25] + fneq[26]);
	P += 2.f*(cxcz * cxcz);

	const float cxcy = (fneq[11] + fneq[12] - fneq[15] - fneq[16] - fneq[19] - fneq[20] + fneq[21] + fneq[22] 
				- fneq[23] - fneq[24] + fneq[25] + fneq[26]);
	P += 2.f*(cxcy * cxcy);

	P = sqrt(P);

	const float CLES_term = 18.f * SmagorinskyConstant * (1.f/rho);
	const float tauLES = 0.5 * tau + 0.5 * sqrt(tau * tau + CLES_term * P);

	omegaLES = 1 / tauLES;
}

__cuda_callable__ void getOuterNormal(const size_t& iCell, const size_t& jCell, const size_t& kCell, int& outerNormalX, int& outerNormalY, int& outerNormalZ, InfoStruct &Info)
{
    outerNormalX = 0;
    outerNormalY = 0;
    outerNormalZ = 0;
    if 			( iCell == 0 ) 						outerNormalX = -1;
    else if 	( iCell == Info.cellCountX - 1 ) 	outerNormalX = 1;
    if 			( jCell == 0 ) 						outerNormalY = -1;
    else if 	( jCell == Info.cellCountY - 1) 	outerNormalY = 1;
    if 			( kCell == 0 ) 						outerNormalZ = -1;
    else if 	( kCell == Info.cellCountZ - 1 ) 	outerNormalZ = 1;
}

__host__ __device__ void getCellIJK(const size_t& cell, size_t& iCell, size_t& jCell, size_t& kCell, InfoStruct &Info)
{
    const size_t xy = Info.cellCountX * Info.cellCountY;
    kCell = cell / xy;
    size_t remainder = cell % xy;
    jCell = remainder / Info.cellCountX;
    iCell = remainder % Info.cellCountX;
}

/*
__host__ __device__ size_t convertIndex(size_t i, size_t j, size_t k, CellCountStruct &cellCount)
{
	size_t cell = k * (cellCount.nx * cellCount.ny) + j * cellCount.nx + i;
	return cell;
}
*/

/*
__host__ __device__ void convertToPhysicalUnits(const float &rho, float &p, float &ux, float &uy, float &uz)
{
	ux = ux * (res/1000.f) / dtPhys;
	uy = uy * (res/1000.f) / dtPhys;
	uz = uz * (res/1000.f) / dtPhys;
	p = (rho - 1.f) * rhoNominalPhys * soundspeedPhys * soundspeedPhys;
}
*/
