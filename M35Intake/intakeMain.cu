constexpr float resGlobal = 2.f; 														// mm
constexpr int gridLevelCount = 3;
constexpr int iterationCount = 20000;
constexpr int iterationChunk = 100;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number
constexpr float rhoOutlet = 1.0f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 20.f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"

#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/applyMirror.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

#include "../STLFunctions.h"
std::string STLPath = "M40IntakeNacaSTL.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	
	if ( Info.gridID == 0 ) // if zero, we are on the coarsest grid
	{
		if ( kCell == 0 || jCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.givenUxUyUz = 1;
		else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
		else Marker.fluid = 1;
	}	
	else // we are on a finer grid
	{
		if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.ghost = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.ghost = 1;
		else if ( kCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
		else Marker.fluid = 1;
	}
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	ux = 0.f;
	uy = 0.f;
	uz = uzInlet;
	const float xPhys = iCell * Info.res + Info.ox;
	const float yPhys = jCell * Info.res + Info.oy;
	const float r2 = xPhys * xPhys + yPhys * yPhys;
	if ( r2 < 18.f * 18.f ) rho = rhoOutlet; // intake
	else rho = 1.f; // lake
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const InfoStruct &Info )
{
	rho = 1.f;
	ux = 0.f;
	uy = 0.f;
	uz = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../exportSectionCutPlot.h"
#include "../fillEquilibrium.h"
#include "../gridRefinementFunctions.h"

void updateGrid( std::vector<GridStruct>& grids, int level ) 
{
    applyStreaming(grids[level]);
    applyLocalCellUpdateBB(grids[level]);
    if (level < gridLevelCount - 1) 
    {
        for (int i = 0; i < 2; i++) updateGrid(grids, level + 1);
        applyCoarseFineGridCommunication(grids[level], grids[level + 1]);
    }
}

int main(int argc, char **argv)
{
	STLStructCPU STLCPU;
	readSTL( STLCPU, STLPath );
	STLStruct STL( STLCPU );
	checkSTLEdges( STL );
	
	std::vector<GridStruct> grids(gridLevelCount);
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	std::cout << "Sizing domain around the STL" << std::endl;
	grids[0].Info.cellCountX = static_cast<int>( std::ceil(( STLCPU.xmax - STLCPU.xmin - 1e-9 ) / grids[0].Info.res ));
	grids[0].Info.cellCountY = static_cast<int>( std::ceil(( STLCPU.ymax - STLCPU.ymin - 1e-9 ) / grids[0].Info.res ));
	grids[0].Info.cellCountZ = static_cast<int>( std::ceil(( STLCPU.zmax - STLCPU.zmin - 1e-9 ) / grids[0].Info.res ));
	grids[0].Info.ox = + STLCPU.xmin + ( 0.5f * ( ( STLCPU.xmax - STLCPU.xmin ) - grids[0].Info.res * ( grids[0].Info.cellCountX-1 ) ) );
	grids[0].Info.oy = + STLCPU.ymin + ( 0.5f * ( ( STLCPU.ymax - STLCPU.ymin ) - grids[0].Info.res * ( grids[0].Info.cellCountY-1 ) ) );
	grids[0].Info.oz = + STLCPU.zmin + ( 0.5f * ( ( STLCPU.zmax - STLCPU.zmin ) - grids[0].Info.res * ( grids[0].Info.cellCountZ-1 ) ) );
	grids[0].Info.cellCountY = grids[0].Info.cellCountY + 1; // adding one more "wall" layer on top
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	grids[0].fArray.setSizes( 27, grids[0].Info.cellCount );
	grids[0].shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( grids[0] );
	grids[0].bouncebackMarkerArray = BoolArrayType( grids[0].Info.cellCount, 0 );
	const bool insideMarkerValue = 0;
	applyMarkersInsideSTL( grids[0].bouncebackMarkerArray, STL, insideMarkerValue, grids[0].Info );
	std::cout << "Cell count on grid " << 0 << ": " << grids[0].Info.cellCount << std::endl;
	
	for ( int level = 1; level < gridLevelCount; level++ )
	{
		grids[level].Info.gridID = grids[level-1].Info.gridID + 1;
		grids[level].Info.res = grids[level-1].Info.res * 0.5f;
		grids[level].Info.dtPhys = grids[level-1].Info.dtPhys * 0.5f;
		grids[level].Info.nu = (grids[level].Info.dtPhys * nuPhys) / ((grids[level].Info.res/1000.f) * (grids[level].Info.res/1000.f));
		
		float progress = (float)level / (float)(gridLevelCount-1);
		progress = std::pow( progress, 0.2f );
		const float xStart = (1 - progress) * grids[0].Info.ox + progress * (-20.f);
		const float xEnd = (1 - progress) * (-grids[0].Info.ox) + progress * 20.f;
		const float yStart = (1 - progress) * grids[0].Info.oy + progress * (-40.f);
		grids[level-1].Info.iSubgridStart = (int)((xStart - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.iSubgridEnd = (int)((xEnd - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.jSubgridStart = (int)((yStart - grids[level-1].Info.oy) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.jSubgridEnd = grids[level-1].Info.cellCountY-1;
		grids[level-1].Info.kSubgridStart = 0;
		grids[level-1].Info.kSubgridEnd = grids[level-1].Info.cellCountZ-1;
		
		grids[level-1].Info.iSubgridStart = std::max({0, grids[level-1].Info.iSubgridStart});
		grids[level-1].Info.iSubgridEnd = std::min({grids[level-1].Info.cellCountX-1, grids[level-1].Info.iSubgridEnd});
		grids[level-1].Info.jSubgridStart = std::max({0, grids[level-1].Info.jSubgridStart});
		grids[level-1].Info.jSubgridEnd = std::min({grids[level-1].Info.cellCountY-1, grids[level-1].Info.jSubgridEnd});
		grids[level-1].Info.kSubgridStart = std::max({0, grids[level-1].Info.kSubgridStart});
		grids[level-1].Info.kSubgridEnd = std::min({grids[level-1].Info.cellCountZ-1, grids[level-1].Info.kSubgridEnd});
		
		grids[level].Info.ox = grids[level-1].Info.ox + grids[level-1].Info.iSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		grids[level].Info.oy = grids[level-1].Info.oy + grids[level-1].Info.jSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		grids[level].Info.oz = grids[level-1].Info.oz + grids[level-1].Info.kSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		
		grids[level].Info.cellCountX = (grids[level-1].Info.iSubgridEnd - grids[level-1].Info.iSubgridStart) * 2;
		grids[level].Info.cellCountY = (grids[level-1].Info.jSubgridEnd - grids[level-1].Info.jSubgridStart) * 2;
		grids[level].Info.cellCountZ = (grids[level-1].Info.kSubgridEnd - grids[level-1].Info.kSubgridStart) * 2;
		grids[level].Info.cellCount = grids[level].Info.cellCountX * grids[level].Info.cellCountY * grids[level].Info.cellCountZ;
		
		grids[level].fArray.setSizes( 27, grids[level].Info.cellCount );	
		grids[level].shifter = IntArrayType( 27, 0 );
		fillEquilibriumFromFunction( grids[level] );
		
		grids[level].bouncebackMarkerArray = BoolArrayType( grids[level].Info.cellCount, 0 );
		applyMarkersInsideSTL( grids[level].bouncebackMarkerArray, STL, insideMarkerValue, grids[level].Info );
		std::cout << "Cell count on grid " << level << ": " << grids[level].Info.cellCount << std::endl;
	}
	
	int cellCountTotal = 0;
	long int cellUpdatesPerIteration = 0;
	for ( int level = 0; level < gridLevelCount; level++ )
	{
		const int cellCountLevel = grids[level].Info.cellCount;
		cellCountTotal += cellCountLevel; 
		cellUpdatesPerIteration += cellCountLevel * std::pow(2, level);
		cellUpdatesPerIteration -= (grids[level].Info.iSubgridEnd - grids[level].Info.iSubgridStart - 2) 
								* (grids[level].Info.jSubgridEnd - grids[level].Info.jSubgridStart - 2) 
								* (grids[level].Info.kSubgridEnd - grids[level].Info.kSubgridStart - 2)
								* std::pow(2, level);
	}
	std::cout << "Cell count total: " << cellCountTotal << std::endl;
	std::cout << "Cell updates per iteration: " << cellUpdatesPerIteration << std::endl;
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			for ( int level = gridLevelCount-2; level >= 0; level-- )
			{
				fillCoarseGridFromFine( grids[level], grids[level+1] );
			}
			
			for ( int level = 0; level < gridLevelCount; level++ )
			{
				const int iCut = grids[level].Info.cellCountX / 2;
				exportSectionCutPlotZY( grids[level], iCut, iteration + level );
			}
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
