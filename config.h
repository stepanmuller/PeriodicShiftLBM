#pragma once

#include <cmath>

// input physics
constexpr float nuPhys = 1.4607e-5;									// m2/s air
constexpr float rhoNominalPhys = 1.225f;								// kg/m3 air
constexpr float res = 0.002f; 										// m
constexpr float RePhys = 6e6; // NOT WORKING CORRECTLY				// 1
constexpr float uzInlet = 0.1f; 										// also works as nominal LBM Mach number
constexpr float SmagorinskyConstant = 0.1f; 							// set to zero to turn off LES

// input domain
size_t cellCountX = 50;
size_t cellCountY = 1000;
size_t cellCountZ = 2000;
constexpr int iterationCount = 2000;

// box dimensions
constexpr size_t boxStartJ = 440;
constexpr size_t boxEndJ = 560;
constexpr size_t boxStartK = 340;
constexpr size_t boxEndK = 460;

// calculated from input
size_t cellCount = cellCountX * cellCountY * cellCountZ;
constexpr float uzInletPhys = nuPhys * RePhys; 						// m/s
constexpr float dtPhys = (uzInlet / uzInletPhys) * res; 			// s
const float soundspeedPhys = (1.f / sqrt(3.f)) * res / dtPhys; 		// m/s
constexpr float nu = (dtPhys * nuPhys) / (res * res);				// LBM nu
constexpr float tau = 3.f * nu + 0.5f;								// LBM tau









