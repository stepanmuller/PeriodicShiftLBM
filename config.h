#pragma once

#include <cmath>

// input physics
constexpr float nuPhys = 1.4607e-5;									// m2/s air
constexpr float rhoNominalPhys = 1.225;								// kg/m3 air
constexpr float res = 0.002; 										// m
constexpr float RePhys = 6e6;										// 1
constexpr float uzInlet = 0.1; 										// also works as nominal LBM Mach number
constexpr float SmagorinskyConstant = 0.1; 							// set to zero to turn off LES

// input domain
constexpr size_t cellCountX = 50;
constexpr size_t cellCountY = 800;
constexpr size_t cellCountZ = 1600;
constexpr int iterationCount = 20000;

// calculated from input
constexpr size_t cellCount = cellCountX * cellCountY * cellCountZ;
constexpr float uzInletPhys = nuPhys * RePhys; 						// m/s given that chord length is 1 m
constexpr float dtPhys = (uzInlet / uzInletPhys) * res; 			// s
const float soundspeedPhys = (1.f / sqrt(3.f)) * res / dtPhys; 		// m/s
constexpr float nu = (dtPhys * nuPhys) / (res * res);				// LBM nu
constexpr float tau = 3.f * nu + 0.5f;								// LBM tau









