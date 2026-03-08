constexpr int caseID = 3;
constexpr float angleOfAttack = 10.f;													// deg

constexpr float chordLengthPhys = 1000.f;												// mm
constexpr float resGlobal = 1.4f; 														// mm
constexpr float uxInlet = 0.05f; 														// also works as nominal LBM Mach number
constexpr float reynoldsNumber = 2000000.f;
constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float nuPhys = 1.50e-5;														// m2/s air
constexpr float uxInletPhys = nuPhys * reynoldsNumber / (chordLengthPhys / 1000.f); 	// m/s
constexpr float rhoNominalPhys = 1.225f;												// kg/m3 air
constexpr float dtPhysGlobal = (uxInlet / uxInletPhys) * (resGlobal/1000); 				// s
constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float domainSizePhys = 10.f * chordLengthPhys;								// mm

const int cellCountX = static_cast<int>(std::ceil(domainSizePhys / resGlobal));
const int cellCountY = cellCountX;
const int cellCountZ = 1;

constexpr int iterationCount = 500000;
constexpr int iterationChunk = 10000;
constexpr int gridLevelCount = 3;

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
std::string STLPathNACA = "NACA0012.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
   	if ( Marker.bounceback ) return;
   	if ( Info.gridID != 0 ) // if not zero, we are on a finer grid
	{
		if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.ghost = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.ghost = 1;
		else if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.periodicZ = 1;
		else Marker.fluid = 1;
	}	
	else // we are on the coarse grid
	{
		Marker.periodicZ = 1;
		if ( iCell == 0 ) Marker.givenUxUyUz = 1;
		else if ( jCell == 0 || jCell == Info.cellCountY-1 ) Marker.givenUxUyUz = 1;
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
	if ( Info.gridID != 0 ) return SmagorinskyConstantGlobal;
	const float x = iCell * Info.res;
	if ( x > domainSizePhys * 0.8f ) return 100.f;
	else return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	rho = 1.f;
	if ( !Marker.bounceback ) ux = uxInlet;
	else ux = 0.f;
	uy = 0.f;
	uz = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../plotter/exportSectionCutPlot.h"
#include "../fillEquilibrium.h"
#include "../gridRefinementFunctions.h"

float getLift( GridStruct &Grid )
{
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	InfoStruct Info = Grid.Info;
	
	const int start = 0;
	const int end = Info.cellCount;
	
	auto fetch = [ = ] __cuda_callable__( const int cell )
	{
		
		MarkerStruct Marker;
		Marker.bounceback = bouncebackMarkerArrayView( cell );
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
		const float gy = rho * (uy -  uyPrev);
		return gy;
	};
	auto reduction = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float gySum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetch, reduction, 0.f );
	float gx, gz = 0;
	convertToPhysicalForce( gx, gySum, gz, Info );
	return gySum;
}

void exportHistoryData( const std::vector<float>& historyVector, const int &currentIteration, int fileNumber ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyVector.data(), sizeof(float), count, fp);
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
	STLStructCPU STLCPUNACA;
	readSTL( STLCPUNACA, STLPathNACA );
	STLStruct STLNACA( STLCPUNACA );
	checkSTLEdges( STLNACA );
	
	float radians = - angleOfAttack * (M_PI / 180.0f);
	rotateSTLAlongZ( STLNACA, radians );
	
	std::vector<GridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = cellCountX;
	grids[0].Info.cellCountY = cellCountY;
	grids[0].Info.cellCountZ = cellCountZ;
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	grids[0].Info.ox = - domainSizePhys * 0.3f;
	grids[0].Info.oy = - domainSizePhys * 0.5f;
	grids[0].fArray.setSizes( 27, grids[0].Info.cellCount );
	grids[0].shifter = IntArrayType( 27, 0 );
	
	grids[0].bouncebackMarkerArray = BoolArrayType( grids[0].Info.cellCount, 0 );
	const bool insideMarkerValue = 1;
	applyMarkersInsideSTL( grids[0].bouncebackMarkerArray, STLNACA, insideMarkerValue, grids[0].Info );
	
	fillEquilibriumFromFunction( grids[0] );
	
	std::cout << "Cell count on grid " << 0 << ": " << grids[0].Info.cellCount << std::endl;
	
	for ( int level = 1; level < gridLevelCount; level++ )
	{
		grids[level].Info.gridID = grids[level-1].Info.gridID + 1;
		grids[level].Info.res = grids[level-1].Info.res * 0.5f;
		grids[level].Info.dtPhys = grids[level-1].Info.dtPhys * 0.5f;
		grids[level].Info.nu = (grids[level].Info.dtPhys * nuPhys) / ((grids[level].Info.res/1000.f) * (grids[level].Info.res/1000.f));
		
		float progress = (float)level / (float)(gridLevelCount-1);
		progress = std::pow( progress, 0.2f );
		
		const float xStart = (1 - progress) * grids[0].Info.ox - progress * 0.6f * chordLengthPhys;
		const float xEnd = (1 - progress) * domainSizePhys + progress * 0.6f * chordLengthPhys;
		grids[level-1].Info.iSubgridStart = (int)((xStart - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.iSubgridEnd = (int)((xEnd - grids[level-1].Info.ox) / grids[level-1].Info.res + 0.5f) + 1;
		
		const float yStart = (1 - progress) * grids[0].Info.oy - progress * 0.3f * chordLengthPhys;
		const float yEnd = - yStart;
		grids[level-1].Info.jSubgridStart = (int)((yStart - grids[level-1].Info.oy) / grids[level-1].Info.res + 0.5f);
		grids[level-1].Info.jSubgridEnd = (int)((yEnd - grids[level-1].Info.oy) / grids[level-1].Info.res + 0.5f) + 1;
		
		grids[level-1].Info.kSubgridStart = 0;
		grids[level-1].Info.kSubgridEnd = grids[level-1].Info.cellCountZ;
		
		grids[level].Info.ox = grids[level-1].Info.ox + grids[level-1].Info.iSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		grids[level].Info.oy = grids[level-1].Info.oy + grids[level-1].Info.jSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		grids[level].Info.oz = grids[level-1].Info.oz + grids[level-1].Info.kSubgridStart * grids[level-1].Info.res - grids[level].Info.res * 0.5f;
		
		grids[level].Info.cellCountX = (grids[level-1].Info.iSubgridEnd - grids[level-1].Info.iSubgridStart) * 2;
		grids[level].Info.cellCountY = (grids[level-1].Info.jSubgridEnd - grids[level-1].Info.jSubgridStart) * 2;
		grids[level].Info.cellCountZ = (grids[level-1].Info.kSubgridEnd - grids[level-1].Info.kSubgridStart) * 2;
		grids[level].Info.cellCount = grids[level].Info.cellCountX * grids[level].Info.cellCountY * grids[level].Info.cellCountZ;
		
		const int cellCountLevel = grids[level].Info.cellCount;
		std::cout << "Cell count on grid " << level << ": " << cellCountLevel << std::endl;
		
		grids[level].fArray.setSizes( 27, grids[level].Info.cellCount );	
		grids[level].shifter = IntArrayType( 27, 0 );
		
		grids[level].bouncebackMarkerArray = BoolArrayType( grids[level].Info.cellCount, 0 );
		const bool insideMarkerValue = 1;
		applyMarkersInsideSTL( grids[level].bouncebackMarkerArray, STLNACA, insideMarkerValue, grids[level].Info );
		
		fillEquilibriumFromFunction( grids[level] );
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
	
	std::vector<float> historyVector( iterationCount, 0.f );
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		
		const float lift = getLift( grids[gridLevelCount-1] );
		const float liftCoefficient = - lift / (0.5 * rhoNominalPhys * uxInletPhys * uxInletPhys * (chordLengthPhys / 1000) * (resGlobal / 1000) * grids[0].Info.cellCountZ);
		
		historyVector[iteration] = liftCoefficient;
		
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
				const int kCut = grids[level].Info.cellCountZ / 2;
				exportSectionCutPlotXY( grids[level], kCut, iteration + level );
				system("python3 ../plotter/plotter.py");
			}
			
			exportHistoryData( historyVector, iteration, caseID );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
