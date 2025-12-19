#pragma once

#include <cmath>

// input physics
constexpr float nuPhys = 1.4607e-5;									// m2/s air
constexpr float rhoNominalPhys = 1.225;								// kg/m3 air
constexpr float res = 0.002; 										// m
constexpr float RePhys = 6e6; // NOT WORKING CORRECTLY				// 1
constexpr float uzInlet = 0.1; 										// also works as nominal LBM Mach number
constexpr float SmagorinskyConstant = 0.1; 							// set to zero to turn off LES

// input domain
constexpr size_t cellCountX = 100;
constexpr size_t cellCountY = 1000;
constexpr size_t cellCountZ = 1000;
constexpr int iterationCount = 10000;

// box dimensions
// constexpr size_t boxStartJ = 350;
// constexpr size_t boxEndJ = 450;
// constexpr size_t boxStartK = 250;
// constexpr size_t boxEndK = 350;

/*
// SMALLER VERSION
// input domain
constexpr size_t cellCountX = 25;
constexpr size_t cellCountY = 400;
constexpr size_t cellCountZ = 800;
constexpr int iterationCount = 10000;

// box dimensions
constexpr size_t boxStartJ = 175;
constexpr size_t boxEndJ = 225;
constexpr size_t boxStartK = 125;
constexpr size_t boxEndK = 175;
*/

// calculated from input
constexpr size_t cellCount = cellCountX * cellCountY * cellCountZ;
constexpr float uzInletPhys = nuPhys * RePhys; 						// m/s
constexpr float dtPhys = (uzInlet / uzInletPhys) * res; 			// s
const float soundspeedPhys = (1.f / sqrt(3.f)) * res / dtPhys; 		// m/s
constexpr float nu = (dtPhys * nuPhys) / (res * res);				// LBM nu
constexpr float tau = 3.f * nu + 0.5f;								// LBM tau









