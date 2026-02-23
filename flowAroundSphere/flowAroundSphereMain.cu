constexpr float sphereDiameterPhys = 2000.f;											// mm
constexpr float resGlobal = 100.f; 														// mm
constexpr float uxInlet = 0.07f; 														// also works as nominal LBM Mach number
constexpr float reynoldsNumber = 1000.f;
constexpr float SmagorinskyConstantGlobal = 0.0f; 										// set to zero to turn off LES

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

constexpr int iterationCount = 100000;
constexpr int iterationChunk = 1000;

#include "../types.h"

#include "../cellFunctions.h"
#include "../applyStreaming.h"
#include "../applyCollision.h"

#include "../multigridFunctions.h"

#include "../boundaryConditions/applyBounceback.h"
#include "../boundaryConditions/applyMirror.h"
#include "../boundaryConditions/restoreRho.h"
#include "../boundaryConditions/restoreUxUyUz.h"
#include "../boundaryConditions/restoreRhoUxUyUz.h"
#include "../boundaryConditions/applyMBBC.h"

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									bool& fluidMarker, bool& bouncebackMarker, bool& mirrorMarker, bool& periodicMarker, bool& givenRhoMarker, bool& givenUxUyUzMarker,
									const InfoStruct& Info )
{
    fluidMarker = 0;
	bouncebackMarker = 0;
	mirrorMarker = 0;
	givenRhoMarker = 0;
	givenUxUyUzMarker = 0;
	
	const float xPhys = iCell * Info.res + Info.ox;
	const float yPhys = jCell * Info.res + Info.oy;
	const float zPhys = kCell * Info.res + Info.oz;
	
	const float r2 = (xPhys-sphereXPhys) * (xPhys - sphereXPhys) + (yPhys - sphereYPhys) * (yPhys - sphereYPhys) + (zPhys - sphereZPhys) * (zPhys - sphereZPhys);
	
	if ( Info.gridID != 0 )
	{
		if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) bouncebackMarker = 1;
		else fluidMarker = 1;
		periodicMarker = 0;
	}	
	else
	{
		if ( r2 <= sphereRadiusPhys * sphereRadiusPhys ) bouncebackMarker = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) givenUxUyUzMarker = 1;
		else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) givenUxUyUzMarker = 1;
		else if ( iCell == 0 ) givenUxUyUzMarker = 1;
		else if ( iCell == Info.cellCountX-1 ) givenRhoMarker = 1;
		else fluidMarker = 1;
		periodicMarker = 0;
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
	bool fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker;
	getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker, Info );
	if ( bouncebackMarker ) ux = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../exportSectionCutPlot.h"
#include "../fillEquilibrium.h"

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
		
		bool fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker;
		getMarkers( iCell, jCell, kCell, fluidMarker, bouncebackMarker, mirrorMarker, periodicMarker, givenRhoMarker, givenUxUyUzMarker, Info );
		
		if ( !bouncebackMarker ) return 0.f;
		
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

void exportHistoryData( const std::vector<float>& historyDragCoefficient, const int &currentIteration ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyDragCoefficient.data(), sizeof(float), count, fp);
    fclose(fp);
    system("python3 historyPlotter.py &"); 
}

int main(int argc, char **argv)
{
	// Coarse grid: Grid0
	GridStruct Grid0;
	Grid0.Info.res = resGlobal;
	Grid0.Info.cellCountX = cellCountX;
	Grid0.Info.cellCountY = cellCountY;
	Grid0.Info.cellCountZ = cellCountZ;
	Grid0.Info.cellCount = Grid0.Info.cellCountX * Grid0.Info.cellCountY * Grid0.Info.cellCountZ;
	std::cout << "Cell count grid 0: " << Grid0.Info.cellCount << std::endl;
	Grid0.Info.dtPhys = dtPhysGlobal;
	Grid0.Info.nu = (Grid0.Info.dtPhys * nuPhys) / ((Grid0.Info.res/1000) * (Grid0.Info.res/1000));
	Grid0.fArray.setSizes( 27, Grid0.Info.cellCount );	
	Grid0.shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( Grid0 );

	Grid0.Info.iSubgridStart = (int)((sphereXPhys - 2.f * sphereRadiusPhys) / Grid0.Info.res);
	Grid0.Info.iSubgridEnd = (int)((sphereXPhys + 2.f * sphereRadiusPhys) / Grid0.Info.res) + 2;
	Grid0.Info.jSubgridStart = (int)((sphereYPhys - 2.f * sphereRadiusPhys) / Grid0.Info.res);
	Grid0.Info.jSubgridEnd = (int)((sphereYPhys + 2.f * sphereRadiusPhys) / Grid0.Info.res) + 2;
	Grid0.Info.kSubgridStart = (int)((sphereZPhys - 2.f * sphereRadiusPhys) / Grid0.Info.res);
	Grid0.Info.kSubgridEnd = (int)((sphereZPhys + 2.f * sphereRadiusPhys) / Grid0.Info.res) + 2;
	
	std::cout << "Grid 0 subgrid bounds: " 
				<< Grid0.Info.iSubgridStart << " " << Grid0.Info.iSubgridEnd << " "
				<< Grid0.Info.jSubgridStart << " " << Grid0.Info.jSubgridEnd << " "
				<< Grid0.Info.kSubgridStart << " " << Grid0.Info.kSubgridEnd << " "
				<< std::endl;
	
	// Finer grid: Grid1
	GridStruct Grid1;
	Grid1.Info.gridID = Grid0.Info.gridID + 1;
	Grid1.Info.res = Grid0.Info.res * 0.5f;
	Grid1.Info.ox = Grid0.Info.iSubgridStart * Grid0.Info.res - Grid1.Info.res * 0.5f;
	Grid1.Info.oy = Grid0.Info.jSubgridStart * Grid0.Info.res - Grid1.Info.res * 0.5f;
	Grid1.Info.oz = Grid0.Info.kSubgridStart * Grid0.Info.res - Grid1.Info.res * 0.5f;
	
	std::cout << "Grid 1 origin: " << Grid1.Info.ox << " " << Grid1.Info.oy << " " << Grid1.Info.oz << " " << std::endl;
	
	Grid1.Info.cellCountX = (Grid0.Info.iSubgridEnd - Grid0.Info.iSubgridStart) * 2;
	Grid1.Info.cellCountY = (Grid0.Info.jSubgridEnd - Grid0.Info.jSubgridStart) * 2;
	Grid1.Info.cellCountZ = (Grid0.Info.kSubgridEnd - Grid0.Info.kSubgridStart) * 2;
	Grid1.Info.cellCount = Grid1.Info.cellCountX * Grid1.Info.cellCountY * Grid1.Info.cellCountZ;
	std::cout << "Cell count grid 1: " << Grid1.Info.cellCount << std::endl;
	Grid1.Info.dtPhys = Grid0.Info.dtPhys * 0.5f;
	Grid1.Info.nu = (Grid1.Info.dtPhys * nuPhys) / ((Grid1.Info.res/1000) * (Grid1.Info.res/1000));
	Grid1.fArray.setSizes( 27, Grid1.Info.cellCount );	
	Grid1.shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( Grid1 );
	
	std::cout << "Starting simulation" << std::endl;
	
	std::vector<float> historyDragCoefficient( iterationCount, 0.f );
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		applyCoarseFineGridCommunication( Grid0, Grid1 );
		
		applyStreaming( Grid0 );
		applyLocalCellUpdate( Grid0 );
		
		for ( int subgridIteration = 0; subgridIteration <= 1; subgridIteration++ )
		{
			applyStreaming( Grid1 );
			applyLocalCellUpdate( Grid1 );
		}
		
		const float drag = getSphereDrag( Grid0 );
		const float dragCoefficient = - (8 * drag) / (rhoNominalPhys * uxInletPhys * uxInletPhys * 3.14159f * (sphereDiameterPhys / 1000.f) * (sphereDiameterPhys / 1000.f));
		
		historyDragCoefficient[iteration] = dragCoefficient;
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			std::cout << "Drag " << drag << std::endl;
			std::cout << "Drag coefficient " << dragCoefficient << std::endl;
			
			const float updateCount = ((float)Grid0.Info.cellCount + 2 * (float)Grid1.Info.cellCount) * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			const int kCut0 = Grid0.Info.cellCountZ / 2;
			exportSectionCutPlotXY( Grid0, kCut0, iteration );
			
			const int kCut1 = Grid1.Info.cellCountZ / 2;
			exportSectionCutPlotXY( Grid1, kCut1, iteration + 1 );
			
			exportHistoryData( historyDragCoefficient, iteration );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
