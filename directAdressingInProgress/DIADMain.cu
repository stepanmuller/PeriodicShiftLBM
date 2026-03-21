constexpr float sphereDiameterPhys = 14.f;											// mm
constexpr float resGlobal = 0.05f; 														// mm
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

constexpr float domainSizePhys = 24.f;													// mm
constexpr float sphereXPhys = 0.5f * domainSizePhys;									// mm
constexpr float sphereYPhys = 0.5f * domainSizePhys;									// mm
constexpr float sphereZPhys = 0.5f * domainSizePhys;									// mm

const int cellCountX = static_cast<int>(std::ceil(domainSizePhys / resGlobal));
const int cellCountY = cellCountX;
const int cellCountZ = cellCountX;

constexpr int iterationCount = 20000;
constexpr int iterationChunk = 1000;
constexpr int gridLevelCount = 4;

constexpr int wallRefinementSpan = 3; // how many cells there are in each refinement layer around the wall

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
#include "./DIADFunctions.h"

//#include "../include/applyLocalCellUpdate.h"
#include "../include/plotter/exportSectionCutPlot.h"
//#include "../include/fillEquilibrium.h"
//#include "../include/gridRefinementFunctions.h"

int main(int argc, char **argv)
{
	std::vector<DIADGridStruct> grids(gridLevelCount);
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = cellCountX;
	grids[0].Info.cellCountY = cellCountY;
	grids[0].Info.cellCountZ = cellCountZ;
	grids[0].Info.ox = 0.5f * grids[0].Info.res;
	grids[0].Info.oy = 0.5f * grids[0].Info.res;
	grids[0].Info.oz = 0.5f * grids[0].Info.res;
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	
	std::cout << "Cell count: " << grids[0].Info.cellCount << std::endl;
	
	buildIJKFromInfo( grids[0].IJK, grids[0].Info );
	
	BoolArrayType fluidMarkerArray = BoolArrayType( grids[0].Info.cellCount );
	std::vector<STLStruct> STLs = {};
	markWhereFinestFluidIs( fluidMarkerArray, grids[0].IJK, STLs, grids[0].Info );
	
	grids[0].bouncebackMarkerArray = fluidMarkerArray;
	const int kCell = grids[0].Info.cellCountZ / 2;
	exportSectionCutPlotXY( grids[0], kCell, 0 );
	system("python3 ../include/plotter/plotter.py");
	
	DIADNeighboursStruct Neighbours;
	
	getDIADNeighbours( Neighbours, grids[0].IJK );
	
	/*
	for ( int cell = 0; cell < grids[0].Info.cellCount; cell++ )
	{
		int iPlus = Neighbours.iPlusArray.getElement(cell);
		int jPlus = Neighbours.jPlusArray.getElement(cell);
		int kPlus = Neighbours.kPlusArray.getElement(cell);
		int iMinus = Neighbours.iMinusArray.getElement(cell);
		int jMinus = Neighbours.jMinusArray.getElement(cell);
		int kMinus = Neighbours.kMinusArray.getElement(cell);		
		
		std::cout << "cell: " << cell <<
		" i: " << grids[0].IJK.iArray.getElement(cell) <<
		" j: " << grids[0].IJK.jArray.getElement(cell) <<
		" k: " << grids[0].IJK.kArray.getElement(cell) << std::endl;
		std::cout << "neighbours: " <<
		" iPlus: " << iPlus <<
		" jPlus: " << jPlus <<
		" kPlus: " << kPlus <<
		" iMinus: " << iMinus <<
		" jMinus: " << jMinus  <<
		" kMinus: " << kMinus  << std::endl;

		iPlus = std::max({0, iPlus});
		jPlus = std::max({0, jPlus});
		kPlus = std::max({0, kPlus});
		iMinus = std::max({0, iMinus});
		jMinus = std::max({0, jMinus});
		kMinus = std::max({0, kMinus});
		std::cout << "neighbour ijk: " << 
		" iPlus i: " << grids[0].IJK.iArray.getElement(iPlus) <<
		" iPlus j: " << grids[0].IJK.jArray.getElement(iPlus) <<
		" iPlus k: " << grids[0].IJK.kArray.getElement(iPlus) <<
		" jPlus i: " << grids[0].IJK.iArray.getElement(jPlus) <<
		" jPlus j: " << grids[0].IJK.jArray.getElement(jPlus) <<
		" jPlus k: " << grids[0].IJK.kArray.getElement(jPlus) <<
		" kPlus i: " << grids[0].IJK.iArray.getElement(kPlus) <<
		" kPlus j: " << grids[0].IJK.jArray.getElement(kPlus) <<
		" kPlus k: " << grids[0].IJK.kArray.getElement(kPlus) <<
		" iMinus i: " << grids[0].IJK.iArray.getElement(iMinus) <<
		" iMinus j: " << grids[0].IJK.jArray.getElement(iMinus) <<
		" iMinus k: " << grids[0].IJK.kArray.getElement(iMinus) <<
		" jMinus i: " << grids[0].IJK.iArray.getElement(jMinus) <<
		" jMinus j: " << grids[0].IJK.jArray.getElement(jMinus) <<
		" jMinus k: " << grids[0].IJK.kArray.getElement(jMinus) <<
		" kMinus i: " << grids[0].IJK.iArray.getElement(kMinus) <<
		" kMinus j: " << grids[0].IJK.jArray.getElement(kMinus) <<
		" kMinus k: " << grids[0].IJK.kArray.getElement(kMinus) <<
		std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

	}
	*/
	
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
