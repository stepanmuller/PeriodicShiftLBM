#include "../config.h"

// Calculate float omegaLES from fneq0, fneq1 ... and rho
// Floats fneq0, fneq1 ... and rho must be defined above

// id: 		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
// cx: 		{ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1 };
// cy: 		{ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1 };
// cz: 		{ 0, 0, 0,-1, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1, 1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1 };

// cx * cx: { 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 };
// cy * cy: { 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
// cz * cz: { 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

// cy * cz: { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1 };
// cx * cz: { 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1,-1,-1, 1, 1 };
// cx * cy: { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,-1,-1, 0, 0,-1,-1, 1, 1,-1,-1, 1, 1 };

float P = 0.f;

const float cxcx = (fneq1 + fneq2 + fneq7 + fneq8 + fneq9 + fneq10 + fneq11 + fneq12 + fneq15 
			+ fneq16 + fneq19 + fneq20 + fneq21 + fneq22 + fneq23 + fneq24 + fneq25 + fneq26);
P += (cxcx * cxcx);

const float cycy = (fneq5 + fneq6 + fneq11 + fneq12 + fneq13 + fneq14 + fneq15 + fneq16 + fneq17
			+ fneq18 + fneq19 + fneq20 + fneq21 + fneq22 + fneq23 + fneq24 + fneq25 + fneq26);
P += (cycy * cycy);

const float czcz = (fneq3 + fneq4 + fneq7 + fneq8 + fneq9 + fneq10 + fneq13 + fneq14 + fneq17
			+ fneq18 + fneq19 + fneq20 + fneq21 + fneq22 + fneq23 + fneq24 + fneq25 + fneq26);
P += (czcz * czcz);

const float cycz = (-fneq13 - fneq14 + fneq17 + fneq18 - fneq19 - fneq20 - fneq21 - fneq22
			+ fneq23 + fneq24 + fneq25 + fneq26);
P += 2.f*(cycz * cycz);

const float cxcz = (-fneq7 - fneq8 + fneq9 + fneq10 + fneq19 + fneq20 - fneq21 - fneq22
			- fneq23 - fneq24 + fneq25 + fneq26);
P += 2.f*(cxcz * cxcz);

const float cxcy = (fneq11 + fneq12 - fneq15 - fneq16 - fneq19 - fneq20 + fneq21 + fneq22 
			- fneq23 - fneq24 + fneq25 + fneq26);
P += 2.f*(cxcy * cxcy);

P = sqrt(P);

const float CLES_term = 18.f * SmagorinskyConstant * (1.f/rho);
const float tauLES = 0.5 * tau + 0.5 * sqrt(tau * tau + CLES_term * P);

const float omegaLES = 1 / tauLES;
