constexpr float sphereDiameterPhys = 2000.f;											// mm
constexpr float resGlobal = 200.f; 														// mm
constexpr float uxInlet = 0.07f; 														// also works as nominal LBM Mach number
float reynoldsNumber = 10000.f;
constexpr float SmagorinskyConstantGlobal = 0.0f; 										// set to zero to turn off LES

constexpr float sphereRadiusPhys = 0.5 * sphereDiameterPhys;							// mm
constexpr float uxInletPhys = uxInlet; 													// m/s, physical velocity set to same as LBM velocity
float nuPhys = uxInletPhys * (sphereDiameterPhys / 1000.f) / reynoldsNumber;	// m2/s
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
constexpr int gridLevelCount = 5;

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

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const InfoStruct &Info )
{
	rho = 1.f;
	ux = uxInlet;
	uy = 0.f;
	uz = 0.f;
	MarkerStruct Marker;
	getMarkers( iCell, jCell, kCell, Marker, Info );
	if ( Marker.bounceback ) ux = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../exportSectionCutPlot.h"
#include "../fillEquilibrium.h"
#include "../gridRefinementFunctions.h"

float getSphereDrag( GridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	InfoStruct Info = Grid.Info;
	
	const int iStart = (int)((sphereXPhys - sphereRadiusPhys - Info.ox) / Info.res);
	const int iEnd = (int)((sphereXPhys + sphereRadiusPhys - Info.ox) / Info.res) + 2;
	const int jStart = (int)((sphereYPhys - sphereRadiusPhys - Info.oy) / Info.res);
	const int jEnd = (int)((sphereYPhys + sphereRadiusPhys - Info.oy) / Info.res) + 2;
	const int kStart = (int)((sphereZPhys - sphereRadiusPhys - Info.oz) / Info.res);
	const int kEnd = (int)((sphereZPhys + sphereRadiusPhys - Info.oz) / Info.res) + 2;
	
	const int start = 0;
	const int end = (iEnd - iStart) * (jEnd - jStart) * (kEnd - kStart);
	
	auto fetch = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iSpan = iEnd - iStart;
		const int jSpan = jEnd - jStart;
		const int kRelative = singleIndex / (iSpan * jSpan);
		const int remainder = singleIndex % (iSpan * jSpan);
		const int jRelative = remainder / iSpan;
		const int iRelative = remainder % iSpan;
		
		const int iCell = iRelative + iStart;
		const int jCell = jRelative + jStart;
		const int kCell = kRelative + kStart;
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		MarkerStruct Marker;
		getMarkers( iCell, jCell, kCell, Marker, Info );
		
		if ( !Marker.bounceback ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz( rho, ux, uy, uz, f );
		applyBounceback( f );
		float rhoPrev, uxPrev, uyPrev, uzPrev;
		getRhoUxUyUz( rhoPrev, uxPrev, uyPrev, uzPrev, f );
		const float gx = rho * (ux -  uxPrev);
		return gx;
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float gxSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0.f );
	float gy, gz = 0;
	convertToPhysicalForce( gxSum, gy, gz, Info );
	return gxSum;
}

void exportHistoryData( const std::vector<float>& historyDragCoefficient, const int &currentIteration, int fileNumber ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyDragCoefficient.data(), sizeof(float), count, fp);
    fclose(fp);
    // Construct the command string to pass the number as an argument
    std::string cmd = "python3 historyPlotter.py " + std::to_string(fileNumber) + " &";
    system(cmd.c_str());
}

void updateGrid( std::vector<GridStruct>& grids, int level ) 
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
	const float reynoldsArray[16] = { 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000 };
	
	for ( int caseIndex = 0; caseIndex < 16; caseIndex++ )
	{
		reynoldsNumber = reynoldsArray[caseIndex];
		nuPhys = uxInletPhys * (sphereDiameterPhys / 1000.f) / reynoldsNumber;
		
		std::vector<GridStruct> grids(gridLevelCount);
		// Coarse grid: Grid0
		grids[0].Info.res = resGlobal;
		grids[0].Info.dtPhys = dtPhysGlobal;
		grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
		grids[0].Info.cellCountX = cellCountX;
		grids[0].Info.cellCountY = cellCountY;
		grids[0].Info.cellCountZ = cellCountZ;
		grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
		grids[0].fArray.setSizes( 27, grids[0].Info.cellCount );
		grids[0].shifter = IntArrayType( 27, 0 );
		fillEquilibriumFromFunction( grids[0] );
		
		for ( int level = 1; level < gridLevelCount; level++ )
		{
			grids[level].Info.gridID = grids[level-1].Info.gridID + 1;
			grids[level].Info.res = grids[level-1].Info.res * 0.5f;
			grids[level].Info.dtPhys = grids[level-1].Info.dtPhys * 0.5f;
			grids[level].Info.nu = (grids[level].Info.dtPhys * nuPhys) / ((grids[level].Info.res/1000.f) * (grids[level].Info.res/1000.f));
			
			float progress = (float)level / (float)(gridLevelCount-1);
			progress = std::pow( progress, 0.5f );
			const float xStart = progress * 1.25f * sphereDiameterPhys;
			const float xEnd = (1 - progress) * domainSizePhys + progress * 2.75f * sphereDiameterPhys;
			const float yStart = progress * 4.75f * sphereDiameterPhys;
			
			grids[level-1].Info.iSubgridStart = (int)((xStart - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f);
			grids[level-1].Info.iSubgridEnd = (int)((xEnd - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f);
			grids[level-1].Info.jSubgridStart = (int)((yStart - grids[level-1].Info.oy) / grids[level-1].Info.res + 0.5f);
			grids[level-1].Info.jSubgridEnd = grids[level-1].Info.jSubgridStart + ( grids[level-1].Info.iSubgridEnd - grids[level-1].Info.iSubgridStart );
			grids[level-1].Info.kSubgridStart = grids[level-1].Info.jSubgridStart;
			grids[level-1].Info.kSubgridEnd = grids[level-1].Info.jSubgridEnd;
			
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
		}
		
		int cellCountTotal = 0;
		int cellUpdatesPerIteration = 0;
		for ( int level = 0; level < gridLevelCount; level++ )
		{
			const int cellCountLevel = grids[level].Info.cellCount;
			std::cout << "Cell count on grid " << level << ": " << cellCountLevel << std::endl;
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
		
		std::vector<float> historyDragCoefficient( iterationCount, 0.f );
		
		TNL::Timer lapTimer;
		lapTimer.reset();
		lapTimer.start();
		for (int iteration=0; iteration<=iterationCount; iteration++)
		{
			updateGrid( grids, 0 );
			
			const float drag = getSphereDrag( grids[gridLevelCount-1] );
			const float dragCoefficient = - (8 * drag) / (rhoNominalPhys * uxInletPhys * uxInletPhys * 3.14159f * (sphereDiameterPhys / 1000.f) * (sphereDiameterPhys / 1000.f));
			
			historyDragCoefficient[iteration] = dragCoefficient;
		}
		lapTimer.stop();
		auto lapTime = lapTimer.getRealTime();
		std::cout << "Finished case with Reynolds number " << reynoldsNumber << std::endl;
		
		const float updateCount = (float)cellUpdatesPerIteration * (float)iterationCount;
		const float glups = updateCount / lapTime / 1000000000.f;
		std::cout << "GLUPS: " << glups << std::endl;
		
		for ( int level = gridLevelCount-2; level >= 0; level-- )
		{
			fillCoarseGridFromFine( grids[level], grids[level+1] );
		}
		
		for ( int level = 0; level < gridLevelCount; level++ )
		{
			const int kCut = grids[level].Info.cellCountZ / 2;
			exportSectionCutPlotXY( grids[level], kCut, reynoldsNumber + level );
		}
		
		exportHistoryData( historyDragCoefficient, iterationCount, reynoldsNumber );
		
		lapTimer.reset();
		lapTimer.start();
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
