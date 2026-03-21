constexpr float sphereDiameterPhys = 8.f;												// mm
constexpr float resGlobal = 0.5f; 														// mm
constexpr float uxInlet = 0.07f; 														// also works as nominal LBM Mach number
constexpr float reynoldsNumber = 1000000.f;
constexpr float SmagorinskyConstantGlobal = 0.0f; 										// set to zero to turn off LES

constexpr float sphereRadiusPhys = 0.5 * sphereDiameterPhys;							// mm
constexpr float uxInletPhys = uxInlet; 													// m/s, physical velocity set to same as LBM velocity
constexpr float nuPhys = uxInletPhys * (sphereDiameterPhys / 1000.f) / reynoldsNumber;	// m2/s
constexpr float rhoNominalPhys = 1.225f;												// kg/m3 air
constexpr float dtPhysGlobal = (uxInlet / uxInletPhys) * (resGlobal/1000); 				// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float domainSizePhys = 16.f;													// mm
constexpr float sphereXPhys = 0.f;														// mm
constexpr float sphereYPhys = 0.f;														// mm
constexpr float sphereZPhys = 0.f;														// mm

const int cellCountX = static_cast<int>(std::ceil(domainSizePhys / resGlobal));
const int cellCountY = cellCountX;
const int cellCountZ = cellCountX;

constexpr int iterationCount = 20000;
constexpr int iterationChunk = 1000;
constexpr int gridLevelCount = 7;

constexpr int wallRefinementSpan = 1; // how many cells there are in each refinement layer around the wall

#include "../include/types.h"

#include "../include/cellFunctions.h"
#include "../include/applyStreaming.h"
#include "../include/applyCollision.h"

#include "../include/boundaryConditions/applyBounceback.h"
#include "../include/boundaryConditions/applyMirror.h"
#include "../include/boundaryConditions/restoreRho.h"
#include "../include/boundaryConditions/restoreUxUyUz.h"
#include "../include/boundaryConditions/restoreRhoUxUyUz.h"
#include "../include/boundaryConditions/applyMBBC.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
   	const float xPhys = iCell * Info.res + Info.ox;
	const float yPhys = jCell * Info.res + Info.oy;
	const float zPhys = kCell * Info.res + Info.oz;
	const float r2 = (xPhys-sphereXPhys) * (xPhys - sphereXPhys) + (yPhys - sphereYPhys) * (yPhys - sphereYPhys) + (zPhys - sphereZPhys) * (zPhys - sphereZPhys);
	if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) Marker.bounceback = 1;
}

#include "../include/STLFunctions.h"
std::string STLPathBlade = "blade.STL";

//#include "../include/applyLocalCellUpdate.h"
#include "../include/plotter/exportSectionCutPlot.h"
//#include "../include/fillEquilibrium.h"
//#include "../include/gridRefinementFunctions.h"

#include "./DIADFunctions.h"

int main(int argc, char **argv)
{
	/*
	STLStructCPU STLCPUBlade;
	readSTL( STLCPUBlade, STLPathBlade );
	STLStruct STLBlade( STLCPUBlade );
	checkSTLEdges( STLBlade );
	float radians = 3.14159f * 2.f * (1.f / 4.f);
	rotateSTLAlongX( STLBlade, radians );
	rotateSTLAlongZ( STLBlade, radians );
	STLBlade.oz = -10.f;
	std::vector<STLStruct> STLs = { STLBlade };
	*/
	std::vector<STLStruct> STLs = { };
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = cellCountX;
	grids[0].Info.cellCountY = cellCountY;
	grids[0].Info.cellCountZ = cellCountZ;
	grids[0].Info.ox = - 0.5f * domainSizePhys;
	grids[0].Info.oy = - 0.5f * domainSizePhys;
	grids[0].Info.oz = - 0.5f * domainSizePhys;
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	buildIJKFromInfo( grids[0].IJK, grids[0].Info );
	
	buildDIADGrid( grids, STLs, 0 );
	
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
