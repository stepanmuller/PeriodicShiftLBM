// Cumulant collision
// Martin Geier, Andrea Pasquali, Martin Schönerr 2017
// Parametrization of the cumulant lattice Boltzman method for fourth order accurate diffusion
// not well conditioned version

// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

// NO FORCING VERSION
__host__ __device__ void applyCollision( float (&f)[27], const float &nu, const float &SmagorinskyConstant )
{
	float rho, ux, uy, uz;
	getRhoUxUyUz( rho, ux, uy, uz, f );
	
	//------------------------------------------------------------------------------------
	// -------------------------------CALCULATING LES OMEGA-------------------------------
	//------------------------------------------------------------------------------------

	float feq[27];
	getFeq(rho, ux, uy, uz, feq);
	float fneq[27];
	getFneq(f, feq, fneq);
	float omegaLES;
	getOmegaLES(fneq, rho, nu, SmagorinskyConstant, omegaLES);
	
	//-------------------------- CUMMULANT COLLISION EQUATIONS ---------------------------
	//------------------------------------------------------------------------------------
	//--------------------------- TRANSFORM TO CENTRAL MOMENTS ---------------------------
	//------------------------------------------------------------------------------------
	
	const float ux2 = ux * ux;
	const float uy2 = uy * uy;
	const float uz2 = uz * uz;
	const float rho_inv = 1.f / rho;
	const float n1o3 = 1.f / 3.f;
	const float n2o3 = 2.f / 3.f;
	const float n4o3 = 4.f / 3.f;
	
	//Eq 6
	const float k_mm0 = (f[21] + f[25]) + f[11];
	const float k_mz0 = (f[8] + f[10]) + f[2];
	const float k_mp0 = (f[24] + f[19]) + f[15];
	const float k_zm0 = (f[14] + f[18]) + f[5];
	const float k_zz0 = (f[4] + f[3]) + f[0];
	const float k_zp0 = (f[17] + f[13]) + f[6];
	const float k_pm0 = (f[20] + f[23]) + f[16];
	const float k_pz0 = (f[9] + f[7]) + f[1];
	const float k_pp0 = (f[26] + f[22]) + f[12];
	
	//Eq 7
	const float k_mm1 = (f[21] - f[25]) - uz * k_mm0;
	const float k_mz1 = (f[8] - f[10]) - uz * k_mz0;
	const float k_mp1 = (f[24] - f[19]) - uz * k_mp0;
	const float k_zm1 = (f[14] - f[18]) - uz * k_zm0;
	const float k_zz1 = (f[4] - f[3]) - uz * k_zz0;
	const float k_zp1 = (f[17] - f[13]) - uz * k_zp0;
	const float k_pm1 = (f[20] - f[23]) - uz * k_pm0;
	const float k_pz1 = (f[9] - f[7]) - uz * k_pz0;
	const float k_pp1 = (f[26] - f[22]) - uz * k_pp0;
	
	//Eq 8
	const float k_mm2 = (f[21] + f[25]) - 2.f * uz * (f[21] - f[25]) + uz2 * k_mm0;
	const float k_mz2 = (f[8] + f[10]) - 2.f * uz * (f[8] - f[10]) + uz2 * k_mz0;
	const float k_mp2 = (f[24] + f[19]) - 2.f * uz * (f[24] - f[19]) + uz2 * k_mp0;
	const float k_zm2 = (f[14] + f[18]) - 2.f * uz * (f[14] - f[18]) + uz2 * k_zm0;
	const float k_zz2 = (f[4] + f[3]) - 2.f * uz * (f[4] - f[3]) + uz2 * k_zz0;
	const float k_zp2 = (f[17] + f[13]) - 2.f * uz * (f[17] - f[13]) + uz2 * k_zp0;
	const float k_pm2 = (f[20] + f[23]) - 2.f * uz * (f[20] - f[23]) + uz2 * k_pm0;
	const float k_pz2 = (f[9] + f[7]) - 2.f * uz * (f[9] - f[7]) + uz2 * k_pz0;
	const float k_pp2 = (f[26] + f[22]) - 2.f * uz * (f[26] - f[22]) + uz2 * k_pp0;

	// Eq 9
	const float k_m00 = (k_mp0 + k_mm0) + k_mz0;
	const float k_z00 = (k_zp0 + k_zm0) + k_zz0;
	const float k_p00 = (k_pp0 + k_pm0) + k_pz0;
	const float k_m01 = (k_mp1 + k_mm1) + k_mz1;
	const float k_z01 = (k_zp1 + k_zm1) + k_zz1;
	const float k_p01 = (k_pp1 + k_pm1) + k_pz1;
	const float k_m02 = (k_mp2 + k_mm2) + k_mz2;
	const float k_z02 = (k_zp2 + k_zm2) + k_zz2;
	const float k_p02 = (k_pp2 + k_pm2) + k_pz2;

	// Eq 10
	const float k_m10 = (k_mp0 - k_mm0) - uy * k_m00;
	const float k_z10 = (k_zp0 - k_zm0) - uy * k_z00;
	const float k_p10 = (k_pp0 - k_pm0) - uy * k_p00;
	const float k_m11 = (k_mp1 - k_mm1) - uy * k_m01;
	const float k_z11 = (k_zp1 - k_zm1) - uy * k_z01;
	const float k_p11 = (k_pp1 - k_pm1) - uy * k_p01;
	const float k_m12 = (k_mp2 - k_mm2) - uy * k_m02;
	const float k_z12 = (k_zp2 - k_zm2) - uy * k_z02;
	const float k_p12 = (k_pp2 - k_pm2) - uy * k_p02;

	// Eq 11
	const float k_m20 = (k_mp0 + k_mm0) - 2.f * uy * (k_mp0 - k_mm0) + uy2 * k_m00;
	const float k_z20 = (k_zp0 + k_zm0) - 2.f * uy * (k_zp0 - k_zm0) + uy2 * k_z00;
	const float k_p20 = (k_pp0 + k_pm0) - 2.f * uy * (k_pp0 - k_pm0) + uy2 * k_p00;
	const float k_m21 = (k_mp1 + k_mm1) - 2.f * uy * (k_mp1 - k_mm1) + uy2 * k_m01;
	const float k_z21 = (k_zp1 + k_zm1) - 2.f * uy * (k_zp1 - k_zm1) + uy2 * k_z01;
	const float k_p21 = (k_pp1 + k_pm1) - 2.f * uy * (k_pp1 - k_pm1) + uy2 * k_p01;
	const float k_m22 = (k_mp2 + k_mm2) - 2.f * uy * (k_mp2 - k_mm2) + uy2 * k_m02;
	const float k_z22 = (k_zp2 + k_zm2) - 2.f * uy * (k_zp2 - k_zm2) + uy2 * k_z02;
	const float k_p22 = (k_pp2 + k_pm2) - 2.f * uy * (k_pp2 - k_pm2) + uy2 * k_p02;

	// Eq 12
	const float k_000 = (k_p00 + k_m00) + k_z00;
	const float k_001 = (k_p01 + k_m01) + k_z01;
	const float C_002 = (k_p02 + k_m02) + k_z02;
	const float C_010 = (k_p10 + k_m10) + k_z10;
	const float C_011 = (k_p11 + k_m11) + k_z11;
	const float C_012 = (k_p12 + k_m12) + k_z12;
	const float C_020 = (k_p20 + k_m20) + k_z20;
	const float C_021 = (k_p21 + k_m21) + k_z21;
	const float k_022 = (k_p22 + k_m22) + k_z22;

	// Eq 13
	const float C_100 = (k_p00 - k_m00) - ux * k_000;
	const float C_101 = (k_p01 - k_m01) - ux * k_001;
	const float C_102 = (k_p02 - k_m02) - ux * C_002;
	const float C_110 = (k_p10 - k_m10) - ux * C_010;
	const float C_111 = (k_p11 - k_m11) - ux * C_011;
	const float k_112 = (k_p12 - k_m12) - ux * C_012;
	const float C_120 = (k_p20 - k_m20) - ux * C_020;
	const float k_121 = (k_p21 - k_m21) - ux * C_021;
	const float k_122 = (k_p22 - k_m22) - ux * k_022;

	// Eq 14
	const float C_200 = (k_p00 + k_m00) - 2.f * ux * (k_p00 - k_m00) + ux2 * k_000;
	const float C_201 = (k_p01 + k_m01) - 2.f * ux * (k_p01 - k_m01) + ux2 * k_001;
	const float k_202 = (k_p02 + k_m02) - 2.f * ux * (k_p02 - k_m02) + ux2 * C_002;
	const float C_210 = (k_p10 + k_m10) - 2.f * ux * (k_p10 - k_m10) + ux2 * C_010;
	const float k_211 = (k_p11 + k_m11) - 2.f * ux * (k_p11 - k_m11) + ux2 * C_011;
	const float k_212 = (k_p12 + k_m12) - 2.f * ux * (k_p12 - k_m12) + ux2 * C_012;
	const float k_220 = (k_p20 + k_m20) - 2.f * ux * (k_p20 - k_m20) + ux2 * C_020;
	const float k_221 = (k_p21 + k_m21) - 2.f * ux * (k_p21 - k_m21) + ux2 * C_021;
	const float k_222 = (k_p22 + k_m22) - 2.f * ux * (k_p22 - k_m22) + ux2 * k_022;
	
	//------------------------------------------------------------------------------------
	//------------------------------ CENTRAL MOM. TO CUMULANTS ---------------------------
	//------------------------------------------------------------------------------------
	
	// Eq 51 from Geier 2015
	const float C_211 = k_211 - (C_200 * C_011 + 2.f * C_101 * C_110) * rho_inv;
	const float C_121 = k_121 - (C_020 * C_101 + 2.f * C_110 * C_011) * rho_inv;
	const float C_112 = k_112 - (C_002 * C_110 + 2.f * C_011 * C_101) * rho_inv;

	// Eq 52 from Geier 2015
	const float C_220 = k_220 - (C_020 * C_200 + 2.f * C_110 * C_110) * rho_inv;
	const float C_022 = k_022 - (C_002 * C_020 + 2.f * C_011 * C_011) * rho_inv;
	const float C_202 = k_202 - (C_200 * C_002 + 2.f * C_101 * C_101) * rho_inv;

	// Eq 53 from Geier 2015
	const float C_122 = k_122 - (C_020 * C_102 + C_002 * C_120 + 4.f * C_011 * C_111 + 2.f * (C_110 * C_012 + C_101 * C_021)) * rho_inv;
	const float C_212 = k_212 - (C_002 * C_210 + C_200 * C_012 + 4.f * C_101 * C_111 + 2.f * (C_011 * C_201 + C_110 * C_102)) * rho_inv;
	const float C_221 = k_221 - (C_200 * C_021 + C_020 * C_201 + 4.f * C_110 * C_111 + 2.f * (C_101 * C_120 + C_011 * C_210)) * rho_inv;

	// Eq 54 from Geier 2015
	const float C_222 = k_222
					  - (4.f * C_111 * C_111 + C_200 * k_022 + C_020 * k_202 + C_002 * k_220
						 + 4.f * (C_011 * k_211 + C_101 * k_121 + C_110 * k_112) + 2.f * (C_120 * C_102 + C_210 * C_012 + C_201 * C_021))
							* rho_inv
					  + (16.f * C_110 * C_101 * C_011 + 4.f * (C_101 * C_101 * C_020 + C_011 * C_011 * C_200 + C_110 * C_110 * C_002)
						 + 2.f * C_200 * C_020 * C_002)
							* rho_inv * rho_inv;

	//------------------------------------------------------------------------------------
	// ------------------------------RELAXATION DEFINITIONS-------------------------------
	//------------------------------------------------------------------------------------
	const float omega1 = omegaLES;	// shear viscosity
	const float omega2 = 1.f;  // bulkViscosity

	const float lambda3 = 0.01f;  // Limiter treshold, section 6 Geier 2017
	const float lambda4 = 0.01f;
	const float lambda5 = 0.01f;
	const float omega3 = 8.f * (omega1 - 2.f) * (omega2 * (3.f * omega1 - 1.f) - 5.f * omega1) 							// eq 111
					   / (8.f * (5.f - 2.f * omega1) * omega1 + omega2 * (8.f + omega1 * (9.f * omega1 - 26.f)));
	const float omega120p102 = omega3 + (1.f - omega3) * fabs(C_120 + C_102) / (rho * lambda3 + fabs(C_120 + C_102));  	// limiter eq 116
	const float omega210p012 = omega3 + (1.f - omega3) * fabs(C_210 + C_012) / (rho * lambda3 + fabs(C_210 + C_012));  	// limiter eq 116
	const float omega201p021 = omega3 + (1.f - omega3) * fabs(C_201 + C_021) / (rho * lambda3 + fabs(C_201 + C_021));  	// limiter eq 116
	const float omega4 = 8.f * (omega1 - 2.f) * (omega1 + omega2 * (3.f * omega1 - 7.f)) 								// eq 112
					   / (omega2 * (56.f - 42.f * omega1 + 9.f * omega1 * omega1) - 8.f * omega1);
	const float omega120m102 = omega4 + (1.f - omega4) * fabs(C_120 - C_102) / (rho * lambda4 + fabs(C_120 - C_102));  	// limiter eq 116
	const float omega210m012 = omega4 + (1.f - omega4) * fabs(C_210 - C_012) / (rho * lambda4 + fabs(C_210 - C_012));  	// limiter eq 116
	const float omega201m021 = omega4 + (1.f - omega4) * fabs(C_201 - C_021) / (rho * lambda4 + fabs(C_201 - C_021));  	// limiter eq 116
	const float omega5 = 																								// eq 113
		24.f * (omega1 - 2.f)
		* (4.f * omega1 * omega1 + omega1 * omega2 * (18.f - 13.f * omega1) + omega2 * omega2 * (2.f + omega1 * (6.f * omega1 - 11.f)))
		/ (16.f * omega1 * omega1 * (omega1 - 6.f) - 2.f * omega1 * omega2 * (216.f + 5.f * omega1 * (9.f * omega1 - 46.f))
		   + omega2 * omega2 * (omega1 * (3.f * omega1 - 10.f) * (15.f * omega1 - 28.f) - 48.f));
	const float omega111 = omega5 + (1.f - omega5) * fabs(C_111) / (rho * lambda5 + fabs(C_111));  						// limiter eq 116
	const float omega6 = 1.f;
	const float omega7 = 1.f;
	const float omega8 = 1.f;
	const float omega9 = 1.f;
	const float omega10 = 1.f;
	
	// Eq 114, 115
	const float A = (4.f * omega1 * omega1 + 2.f * omega1 * omega2 * (omega1 - 6.f) + omega2 * omega2 * (omega1 * (10.f - 3.f * omega1) - 4.f))
				  / (omega1 - omega2) / (omega2 * (2.f + 3.f * omega1) - 8.f * omega1);
	const float B =
		(4.f * omega1 * omega2 * (9.f * omega1 - 16.f) - 4.f * omega1 * omega1 - 2.f * omega2 * omega2 * (2.f + 9.f * omega1 * (omega1 - 2.f)))
		* n1o3 / (omega1 - omega2) / (omega2 * (2.f + 3.f * omega1) - 8.f * omega1);

	//------------------------------------------------------------------------------------
	// -------------------------------------COLLISION-------------------------------------
	//------------------------------------------------------------------------------------
	
	const float Cs_110 = (1.f - omega1) * C_110;
	const float Cs_101 = (1.f - omega1) * C_101;
	const float Cs_011 = (1.f - omega1) * C_011;

	// Eq 27-29
	const float Dxu =
		-omega1 * 0.5f * rho_inv * (2.f * C_200 - C_020 - C_002)
		- omega2 * 0.5f * rho_inv * (C_200 + C_020 + C_002 - k_000);
	const float Dyv = Dxu + 1.5f * omega1 * rho_inv * (C_200 - C_020);
	const float Dzw = Dxu + 1.5f * omega1 * rho_inv * (C_200 - C_002);
	
	// Eq 30 - 32
	const float DxvDyu = -3.f * omega1 * rho_inv * C_110;
	const float DxwDzu = -3.f * omega1 * rho_inv * C_101;
	const float DywDzv = -3.f * omega1 * rho_inv * C_011;

	// Eq 33-35
	const float Eq33RHS = (1.f - omega1) * (C_200 - C_020) - 3.f * rho * (1.f - omega1 * 0.5f) * (ux2 * Dxu - uy2 * Dyv);
	const float Eq34RHS = (1.f - omega1) * (C_200 - C_002) - 3.f * rho * (1.f - omega1 * 0.5f) * (ux2 * Dxu - uy2 * Dzw);
	const float Eq35RHS = k_000 * omega2 + (1.f - omega2) * (C_200 + C_020 + C_002)
						- 3.f * rho * (1.f - omega2 * 0.5f) * (ux2 * Dxu + uy2 * Dyv + uy2 * Dzw);

	const float Cs_200 = n1o3 * (Eq33RHS + Eq34RHS + Eq35RHS);
	const float Cs_020 = n1o3 * (-2.f * Eq33RHS + Eq34RHS + Eq35RHS);
	const float Cs_002 = n1o3 * (Eq33RHS - 2.f * Eq34RHS + Eq35RHS);

	// Eq 117-122
	const float Eq117 = (1.f - omega120p102) * (C_120 + C_102);
	const float Eq118 = (1.f - omega210p012) * (C_210 + C_012);
	const float Eq119 = (1.f - omega201p021) * (C_201 + C_021);
	const float Eq120 = (1.f - omega120m102) * (C_120 - C_102);
	const float Eq121 = (1.f - omega210m012) * (C_210 - C_012);
	const float Eq122 = (1.f - omega201m021) * (C_201 - C_021);

	const float Cs_120 = 0.5f * (Eq120 + Eq117);
	const float Cs_102 = 0.5f * (-Eq120 + Eq117);
	const float Cs_210 = 0.5f * (Eq121 + Eq118);
	const float Cs_012 = 0.5f * (-Eq121 + Eq118);
	const float Cs_021 = 0.5f * (-Eq122 + Eq119);
	const float Cs_201 = 0.5f * (Eq122 + Eq119);
	
	// Eq 42
	const float Cs_111 = (1.f - omega111) * C_111;

	// Eq 43-45
	const float Eq43RHS =
		n2o3 * (1.f / omega1 - 0.5f) * omega6 * A * rho * (Dxu - 2.f * Dyv + Dzw) + (1.f - omega6) * (C_220 - 2.f * C_202 + C_022);
	const float Eq44RHS =
		n2o3 * (1.f / omega1 - 0.5f) * omega6 * A * rho * (Dxu + Dyv - 2.f * Dzw) + (1.f - omega6) * (C_220 + C_202 - 2.f * C_022);
	const float Eq45RHS = -n4o3 * (1.f / omega1 - 0.5f) * omega7 * A * rho * (Dxu + Dyv + Dzw) + (1.f - omega7) * (C_220 + C_202 + C_022);

	const float Cs_220 = n1o3 * (Eq43RHS + Eq44RHS + Eq45RHS);
	const float Cs_202 = n1o3 * (-Eq43RHS + Eq45RHS);
	const float Cs_022 = n1o3 * (-Eq44RHS + Eq45RHS);
	// Eq 46-48
	const float Cs_211 = -n1o3 * (1.f / omega1 - 0.5f) * omega8 * B * rho * DywDzv + (1.f - omega8) * C_211;
	const float Cs_121 = -n1o3 * (1.f / omega1 - 0.5f) * omega8 * B * rho * DxwDzu + (1.f - omega8) * C_121;
	const float Cs_112 = -n1o3 * (1.f / omega1 - 0.5f) * omega8 * B * rho * DxvDyu + (1.f - omega8) * C_112;
	// Eqs 49-52
	const float Cs_221 = (1.f - omega9) * C_221;
	const float Cs_212 = (1.f - omega9) * C_212;
	const float Cs_122 = (1.f - omega9) * C_122;
	const float Cs_222 = (1.f - omega10) * C_222;

	//------------------------------------------------------------------------------------
	//------------------------------ CUMULANTS TO CENTRAL MOM. ---------------------------
	//------------------------------------------------------------------------------------
	
	// Eq 81 from Geier 2015
	const float ks_211 = Cs_211 + (Cs_200 * Cs_011 + 2.f * Cs_101 * Cs_110) * rho_inv;
	const float ks_121 = Cs_121 + (Cs_020 * Cs_101 + 2.f * Cs_110 * Cs_011) * rho_inv;
	const float ks_112 = Cs_112 + (Cs_002 * Cs_110 + 2.f * Cs_011 * Cs_101) * rho_inv;

	// Eq 82 from Geier 2015
	const float ks_220 = Cs_220 + (Cs_020 * Cs_200 + 2.f * Cs_110 * Cs_110) * rho_inv;
	const float ks_022 = Cs_022 + (Cs_002 * Cs_020 + 2.f * Cs_011 * Cs_011) * rho_inv;
	const float ks_202 = Cs_202 + (Cs_200 * Cs_002 + 2.f * Cs_101 * Cs_101) * rho_inv;

	// Eq 83 from Geier 2015
	const float ks_122 =
		Cs_122 + (Cs_020 * Cs_102 + Cs_002 * Cs_120 + 4.f * Cs_011 * Cs_111 + 2.f * (Cs_110 * Cs_012 + Cs_101 * Cs_021)) * rho_inv;
	const float ks_212 =
		Cs_212 + (Cs_002 * Cs_210 + Cs_200 * Cs_012 + 4.f * Cs_101 * Cs_111 + 2.f * (Cs_011 * Cs_201 + Cs_110 * Cs_102)) * rho_inv;
	const float ks_221 =
		Cs_221 + (Cs_200 * Cs_021 + Cs_020 * Cs_201 + 4.f * Cs_110 * Cs_111 + 2.f * (Cs_101 * Cs_120 + Cs_011 * Cs_210)) * rho_inv;

	// Eq 84 from Geier 2015
	const float ks_222 =
		Cs_222
		+ (4.f * Cs_111 * Cs_111 + Cs_200 * ks_022 + Cs_020 * ks_202 + Cs_002 * ks_220
		   + 4.f * (Cs_011 * ks_211 + Cs_101 * ks_121 + Cs_110 * ks_112) + 2.f * (Cs_120 * Cs_102 + Cs_210 * Cs_012 + Cs_201 * Cs_021))
			  * rho_inv
		- (16.f * Cs_110 * Cs_101 * Cs_011 + 4.f * (Cs_101 * Cs_101 * Cs_020 + Cs_011 * Cs_011 * Cs_200 + Cs_110 * Cs_110 * Cs_002)
		   + 2.f * Cs_200 * Cs_020 * Cs_002)
			  * rho_inv * rho_inv;

	const float ks_000 = k_000;
	
	// Geier 2017: forcing scheme
	const float Cs_100 = -C_100;
	const float Cs_010 = -C_010;
	const float Cs_001 = -k_001;

	// Eq 88 from Geier 2015
	const float ks_z00 = ks_000 * (1.f - ux2) - 2.f * ux * Cs_100 - Cs_200;
	const float ks_z01 = Cs_001 * (1.f - ux2) - 2.f * ux * Cs_101 - Cs_201;
	const float ks_z02 = Cs_002 * (1.f - ux2) - 2.f * ux * Cs_102 - ks_202;
	const float ks_z10 = Cs_010 * (1.f - ux2) - 2.f * ux * Cs_110 - Cs_210;
	const float ks_z11 = Cs_011 * (1.f - ux2) - 2.f * ux * Cs_111 - ks_211;
	const float ks_z12 = Cs_012 * (1.f - ux2) - 2.f * ux * ks_112 - ks_212;
	const float ks_z20 = Cs_020 * (1.f - ux2) - 2.f * ux * Cs_120 - ks_220;
	const float ks_z21 = Cs_021 * (1.f - ux2) - 2.f * ux * ks_121 - ks_221;
	const float ks_z22 = ks_022 * (1.f - ux2) - 2.f * ux * ks_122 - ks_222;

	// Eq 89 from Geier 2015
	const float ks_m00 = (ks_000 * (ux2 - ux) + Cs_100 * (2.f * ux - 1.f) + Cs_200) * 0.5f;
	const float ks_m01 = (Cs_001 * (ux2 - ux) + Cs_101 * (2.f * ux - 1.f) + Cs_201) * 0.5f;
	const float ks_m02 = (Cs_002 * (ux2 - ux) + Cs_102 * (2.f * ux - 1.f) + ks_202) * 0.5f;
	const float ks_m10 = (Cs_010 * (ux2 - ux) + Cs_110 * (2.f * ux - 1.f) + Cs_210) * 0.5f;
	const float ks_m11 = (Cs_011 * (ux2 - ux) + Cs_111 * (2.f * ux - 1.f) + ks_211) * 0.5f;
	const float ks_m12 = (Cs_012 * (ux2 - ux) + ks_112 * (2.f * ux - 1.f) + ks_212) * 0.5f;
	const float ks_m20 = (Cs_020 * (ux2 - ux) + Cs_120 * (2.f * ux - 1.f) + ks_220) * 0.5f;
	const float ks_m21 = (Cs_021 * (ux2 - ux) + ks_121 * (2.f * ux - 1.f) + ks_221) * 0.5f;
	const float ks_m22 = (ks_022 * (ux2 - ux) + ks_122 * (2.f * ux - 1.f) + ks_222) * 0.5f;

	// Eq 90 from Geier 2015
	const float ks_p00 = (ks_000 * (ux2 + ux) + Cs_100 * (2.f * ux + 1.f) + Cs_200) * 0.5f;
	const float ks_p01 = (Cs_001 * (ux2 + ux) + Cs_101 * (2.f * ux + 1.f) + Cs_201) * 0.5f;
	const float ks_p02 = (Cs_002 * (ux2 + ux) + Cs_102 * (2.f * ux + 1.f) + ks_202) * 0.5f;
	const float ks_p10 = (Cs_010 * (ux2 + ux) + Cs_110 * (2.f * ux + 1.f) + Cs_210) * 0.5f;
	const float ks_p11 = (Cs_011 * (ux2 + ux) + Cs_111 * (2.f * ux + 1.f) + ks_211) * 0.5f;
	const float ks_p12 = (Cs_012 * (ux2 + ux) + ks_112 * (2.f * ux + 1.f) + ks_212) * 0.5f;
	const float ks_p20 = (Cs_020 * (ux2 + ux) + Cs_120 * (2.f * ux + 1.f) + ks_220) * 0.5f;
	const float ks_p21 = (Cs_021 * (ux2 + ux) + ks_121 * (2.f * ux + 1.f) + ks_221) * 0.5f;
	const float ks_p22 = (ks_022 * (ux2 + ux) + ks_122 * (2.f * ux + 1.f) + ks_222) * 0.5f;

	// Eq 91 from Geier 2015
	const float ks_mz0 = ks_m00 * (1.f - uy2) - 2.f * uy * ks_m10 - ks_m20;
	const float ks_mz1 = ks_m01 * (1.f - uy2) - 2.f * uy * ks_m11 - ks_m21;
	const float ks_mz2 = ks_m02 * (1.f - uy2) - 2.f * uy * ks_m12 - ks_m22;
	const float ks_zz0 = ks_z00 * (1.f - uy2) - 2.f * uy * ks_z10 - ks_z20;
	const float ks_zz1 = ks_z01 * (1.f - uy2) - 2.f * uy * ks_z11 - ks_z21;
	const float ks_zz2 = ks_z02 * (1.f - uy2) - 2.f * uy * ks_z12 - ks_z22;
	const float ks_pz0 = ks_p00 * (1.f - uy2) - 2.f * uy * ks_p10 - ks_p20;
	const float ks_pz1 = ks_p01 * (1.f - uy2) - 2.f * uy * ks_p11 - ks_p21;
	const float ks_pz2 = ks_p02 * (1.f - uy2) - 2.f * uy * ks_p12 - ks_p22;

	// Eq 92 from Geier 2015
	const float ks_mm0 = (ks_m00 * (uy2 - uy) + ks_m10 * (2.f * uy - 1.f) + ks_m20) * 0.5f;
	const float ks_mm1 = (ks_m01 * (uy2 - uy) + ks_m11 * (2.f * uy - 1.f) + ks_m21) * 0.5f;
	const float ks_mm2 = (ks_m02 * (uy2 - uy) + ks_m12 * (2.f * uy - 1.f) + ks_m22) * 0.5f;
	const float ks_zm0 = (ks_z00 * (uy2 - uy) + ks_z10 * (2.f * uy - 1.f) + ks_z20) * 0.5f;
	const float ks_zm1 = (ks_z01 * (uy2 - uy) + ks_z11 * (2.f * uy - 1.f) + ks_z21) * 0.5f;
	const float ks_zm2 = (ks_z02 * (uy2 - uy) + ks_z12 * (2.f * uy - 1.f) + ks_z22) * 0.5f;
	const float ks_pm0 = (ks_p00 * (uy2 - uy) + ks_p10 * (2.f * uy - 1.f) + ks_p20) * 0.5f;
	const float ks_pm1 = (ks_p01 * (uy2 - uy) + ks_p11 * (2.f * uy - 1.f) + ks_p21) * 0.5f;
	const float ks_pm2 = (ks_p02 * (uy2 - uy) + ks_p12 * (2.f * uy - 1.f) + ks_p22) * 0.5f;

	// Eq 93 from Geier 2015
	const float ks_mp0 = (ks_m00 * (uy2 + uy) + ks_m10 * (2.f * uy + 1.f) + ks_m20) * 0.5f;
	const float ks_mp1 = (ks_m01 * (uy2 + uy) + ks_m11 * (2.f * uy + 1.f) + ks_m21) * 0.5f;
	const float ks_mp2 = (ks_m02 * (uy2 + uy) + ks_m12 * (2.f * uy + 1.f) + ks_m22) * 0.5f;
	const float ks_zp0 = (ks_z00 * (uy2 + uy) + ks_z10 * (2.f * uy + 1.f) + ks_z20) * 0.5f;
	const float ks_zp1 = (ks_z01 * (uy2 + uy) + ks_z11 * (2.f * uy + 1.f) + ks_z21) * 0.5f;
	const float ks_zp2 = (ks_z02 * (uy2 + uy) + ks_z12 * (2.f * uy + 1.f) + ks_z22) * 0.5f;
	const float ks_pp0 = (ks_p00 * (uy2 + uy) + ks_p10 * (2.f * uy + 1.f) + ks_p20) * 0.5f;
	const float ks_pp1 = (ks_p01 * (uy2 + uy) + ks_p11 * (2.f * uy + 1.f) + ks_p21) * 0.5f;
	const float ks_pp2 = (ks_p02 * (uy2 + uy) + ks_p12 * (2.f * uy + 1.f) + ks_p22) * 0.5f;
	
	//------------------------------------------------------------------------------------
	//----------------------- TRANSFORM TO DISTRIBUTION FUNCTION -------------------------
	//------------------------------------------------------------------------------------
	
	// Eq 94 Geiger 2015
	f[11] = ks_mm0 * (1.f - uz2) - 2.f * uz * ks_mm1 - ks_mm2;
	f[2] = ks_mz0 * (1.f - uz2) - 2.f * uz * ks_mz1 - ks_mz2;
	f[15] = ks_mp0 * (1.f - uz2) - 2.f * uz * ks_mp1 - ks_mp2;
	f[5] = ks_zm0 * (1.f - uz2) - 2.f * uz * ks_zm1 - ks_zm2;
	f[0] = ks_zz0 * (1.f - uz2) - 2.f * uz * ks_zz1 - ks_zz2;
	f[6] = ks_zp0 * (1.f - uz2) - 2.f * uz * ks_zp1 - ks_zp2;
	f[16] = ks_pm0 * (1.f - uz2) - 2.f * uz * ks_pm1 - ks_pm2;
	f[1] = ks_pz0 * (1.f - uz2) - 2.f * uz * ks_pz1 - ks_pz2;
	f[12] = ks_pp0 * (1.f - uz2) - 2.f * uz * ks_pp1 - ks_pp2;

	// Eq 95 Geiger 2015
	f[25] = (ks_mm0 * (uz2 - uz) + ks_mm1 * (2.f * uz - 1.f) + ks_mm2) * 0.5f;
	f[10] = (ks_mz0 * (uz2 - uz) + ks_mz1 * (2.f * uz - 1.f) + ks_mz2) * 0.5f;
	f[19] = (ks_mp0 * (uz2 - uz) + ks_mp1 * (2.f * uz - 1.f) + ks_mp2) * 0.5f;
	f[18] = (ks_zm0 * (uz2 - uz) + ks_zm1 * (2.f * uz - 1.f) + ks_zm2) * 0.5f;
	f[3] = (ks_zz0 * (uz2 - uz) + ks_zz1 * (2.f * uz - 1.f) + ks_zz2) * 0.5f;
	f[13] = (ks_zp0 * (uz2 - uz) + ks_zp1 * (2.f * uz - 1.f) + ks_zp2) * 0.5f;
	f[23] = (ks_pm0 * (uz2 - uz) + ks_pm1 * (2.f * uz - 1.f) + ks_pm2) * 0.5f;
	f[7] = (ks_pz0 * (uz2 - uz) + ks_pz1 * (2.f * uz - 1.f) + ks_pz2) * 0.5f;
	f[22] = (ks_pp0 * (uz2 - uz) + ks_pp1 * (2.f * uz - 1.f) + ks_pp2) * 0.5f;

	// Eq 96 Geiger 2015
	f[21] = (ks_mm0 * (uz2 + uz) + ks_mm1 * (2.f * uz + 1.f) + ks_mm2) * 0.5f;
	f[8] = (ks_mz0 * (uz2 + uz) + ks_mz1 * (2.f * uz + 1.f) + ks_mz2) * 0.5f;
	f[24] = (ks_mp0 * (uz2 + uz) + ks_mp1 * (2.f * uz + 1.f) + ks_mp2) * 0.5f;
	f[14] = (ks_zm0 * (uz2 + uz) + ks_zm1 * (2.f * uz + 1.f) + ks_zm2) * 0.5f;
	f[4] = (ks_zz0 * (uz2 + uz) + ks_zz1 * (2.f * uz + 1.f) + ks_zz2) * 0.5f;
	f[17] = (ks_zp0 * (uz2 + uz) + ks_zp1 * (2.f * uz + 1.f) + ks_zp2) * 0.5f;
	f[20] = (ks_pm0 * (uz2 + uz) + ks_pm1 * (2.f * uz + 1.f) + ks_pm2) * 0.5f;
	f[9] = (ks_pz0 * (uz2 + uz) + ks_pz1 * (2.f * uz + 1.f) + ks_pz2) * 0.5f;
	f[26] = (ks_pp0 * (uz2 + uz) + ks_pp1 * (2.f * uz + 1.f) + ks_pp2) * 0.5f;
}

/*
// FORCING VERSION
__host__ __device__ void applyCollision(
	float &rho, float &ux, float &uy, float &uz, 
	float (&f)[27]
	)
{
	//------------------------------------------------------------------------------------
	//---------------------------- APPLY FORCING - FIRST HALF ----------------------------
	//------------------------------------------------------------------------------------

	const float gx  = 0; //gxArrayView[cell];
	const float gy  = 0; //gyArrayView[cell];
	const float gz  = 0; //gzArrayView[cell];
	ux = ((ux * rho) + gx/2.f) / rho;
	uy = ((uy * rho) + gy/2.f) / rho;
	uz = ((uz * rho) + gz/2.f) / rho;

	//-------------------------- CUMMULANT COLLISION EQUATIONS ---------------------------
	//------------------------------------------------------------------------------------
	//--------------------------- TRANSFORM TO CENTRAL MOMENTS ---------------------------
	//------------------------------------------------------------------------------------

	//Eq Geiger 2015(43)
	//first part of the central moments transformation
	const float k_aa0 = (f[21] + f[25]) + f[11];
	const float k_ab0 = (f[8] + f[10]) + f[2];
	const float k_ac0 = (f[24] + f[19]) + f[15];
	const float k_ba0 = (f[14] + f[18]) + f[5];
	const float k_bb0 = (f[4] + f[3]) + f[0];
	const float k_bc0 = (f[17] + f[13]) + f[6];
	const float k_ca0 = (f[20] + f[23]) + f[16];
	const float k_cb0 = (f[9] + f[7]) + f[1];
	const float k_cc0 = (f[26] + f[22]) + f[12];
(...)
	//Eq  Geiger 2015(96)
	f[21] = (ks_aa0 * (uz2 + uz) + ks_aa1 * (2.f * uz + 1.f) + ks_aa2) * 0.5f;
	f[8] = (ks_ab0 * (uz2 + uz) + ks_ab1 * (2.f * uz + 1.f) + ks_ab2) * 0.5f;
	f[24] = (ks_ac0 * (uz2 + uz) + ks_ac1 * (2.f * uz + 1.f) + ks_ac2) * 0.5f;
	f[14] = (ks_ba0 * (uz2 + uz) + ks_ba1 * (2.f * uz + 1.f) + ks_ba2) * 0.5f;
	f[4] = (ks_bb0 * (uz2 + uz) + ks_bb1 * (2.f * uz + 1.f) + ks_bb2) * 0.5f;
	f[17] = (ks_bc0 * (uz2 + uz) + ks_bc1 * (2.f * uz + 1.f) + ks_bc2) * 0.5f;
	f[20] = (ks_ca0 * (uz2 + uz) + ks_ca1 * (2.f * uz + 1.f) + ks_ca2) * 0.5f;
	f[9] = (ks_cb0 * (uz2 + uz) + ks_cb1 * (2.f * uz + 1.f) + ks_cb2) * 0.5f;
	f[26] = (ks_cc0 * (uz2 + uz) + ks_cc1 * (2.f * uz + 1.f) + ks_cc2) * 0.5f;

	//------------------------------------------------------------------------------------
	//------------------------- FINISH FORCING - SECOND HALF -----------------------------
	//------------------------------------------------------------------------------------

	ux = ((ux * rho) + gx/2.f) / rho;
	uy = ((uy * rho) + gy/2.f) / rho;
	uz = ((uz * rho) + gz/2.f) / rho;
}
*/
