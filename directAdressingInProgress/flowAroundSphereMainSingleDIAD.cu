constexpr float sphereDiameterPhys = 2000.f;											// mm
constexpr float resGlobal = 200.f; 														// mm
constexpr float uxInlet = 0.07f; 														// also works as nominal LBM Mach number
constexpr float reynoldsNumber = 1000000.f;
constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float sphereRadiusPhys = 0.5 * sphereDiameterPhys;							// mm
constexpr float uxInletPhys = uxInlet; 													// m/s, physical velocity set to same as LBM velocity
constexpr float nuPhys = uxInletPhys * (sphereDiameterPhys / 1000.f) / reynoldsNumber;	// m2/s
constexpr float rhoNominalPhys = 1.225f;												// kg/m3 air
constexpr float dtPhysGlobal = (uxInlet / uxInletPhys) * (resGlobal/1000); 				// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float domainSizePhys = 11.f * sphereDiameterPhys;								// mm
constexpr float sphereXPhys = 2.f * sphereDiameterPhys;									// mm
constexpr float sphereYPhys = 5.5f * sphereDiameterPhys;								// mm
constexpr float sphereZPhys = 5.5f * sphereDiameterPhys;								// mm

const int cellCountX = static_cast<int>(std::ceil(domainSizePhys / resGlobal));
const int cellCountY = cellCountX;
const int cellCountZ = cellCountX;

constexpr int iterationCount = 20000;
constexpr int iterationChunk = 1000;
constexpr int gridLevelCount = 5;
constexpr int wallRefinementSpan = 2;

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
	
	if ( Info.gridID != 0 ) // if not zero, we are on a finer grid
	{
		if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) Marker.bounceback = 1;
		else if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.ghost = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.ghost = 1;
		else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.ghost = 1;
		else Marker.fluid = 1;
	}	
	else // we are on the coarse grid
	{
		if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) Marker.bounceback = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.givenUxUyUz = 1;
		else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.givenUxUyUz = 1;
		else if ( iCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( iCell == Info.cellCountX-1 ) Marker.givenRho = 1;
		else Marker.fluid = 1;
	}
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											const InfoStruct& Info )
{
    rho = 1.f;
	ux = uxInlet;
	uy = 0.f;
	uz = 0.f;
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	rho = 1.f;
	if ( Marker.bounceback ) ux = uxInlet;
	else ux = 0.f;
	uy = 0.f;
	uz = 0.f;
}

#include "../include/applyLocalCellUpdate.h"
#include "../include/plotter/exportSectionCutPlot.h"
#include "../include/fillEquilibrium.h"
#include "../include/gridRefinementFunctions.h"

#include "./DIADFunctions.h"

void updateGrid( std::vector<DIADGridStruct>& grids, int level ) 
{
    applyStreaming(grids[level]);
    applyLocalCellUpdate(grids[level]);
    if (level < gridLevelCount - 1) 
    {
        for (int i = 0; i < 2; i++) updateGrid(grids, level + 1);
        applyCoarseFineGridCommunication(grids[level], grids[level + 1]);
    }
}

int main(int argc, char **argv)
{
	std::vector<STLStruct> STLs = { };
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = cellCountX;
	grids[0].Info.cellCountY = cellCountY;
	grids[0].Info.cellCountZ = cellCountZ;
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	buildIJKFromInfo( grids[0].IJK, grids[0].Info );
	
	TNL::Timer buildIJKTimer;
	buildIJKTimer.reset();
	buildIJKTimer.start();
	buildDIADGrids( grids, STLs, 0 );
	buildIJKTimer.stop();
	auto buildIJKTime = buildIJKTimer.getRealTime();
	std::cout << "This took " << buildIJKTime << " s" << std::endl; //7.11 before find IJK rewrite
	
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		grids[level].fArray.setSizes( 27, grids[level].Info.cellCount );
		fillEquilibriumFromFunction( grids[level] );
	}
	
	int cellCountTotal = 0;
	long int cellUpdatesPerIteration = 0;
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		const int cellCountLevel = grids[level].Info.cellCount;
		cellCountTotal += cellCountLevel; 
		cellUpdatesPerIteration += cellCountLevel * std::pow(2, level);
	}
	std::cout << "Cell count total: " << cellCountTotal << std::endl;
	std::cout << "Cell updates per iteration: " << cellUpdatesPerIteration << std::endl;
	
	std::cout << "Starting simulation" << std::endl;
	
	//std::vector<float> historyDragCoefficient( iterationCount, 0.f );
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		
		//const float drag = getSphereDrag( grids[gridLevelCount-1] );
		//const float dragCoefficient = - (8 * drag) / (rhoNominalPhys * uxInletPhys * uxInletPhys * 3.14159f * (sphereDiameterPhys / 1000.f) * (sphereDiameterPhys / 1000.f));
		
		//historyDragCoefficient[iteration] = dragCoefficient;
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			/*
			for ( int level = gridLevelCount-2; level >= 0; level-- )
			{
				fillCoarseGridFromFine( grids[level], grids[level+1] );
			}
			*/
			
			for ( int level = 0; level < gridLevelCount; level++ )
			{
				const int kCut = grids[level].Info.cellCountZ / 2;
				exportSectionCutPlotXY( grids[level], kCut, iteration + level );
				system("python3 ../include/plotter/plotter.py");
			}
			
			//exportHistoryData( historyDragCoefficient, iteration );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
