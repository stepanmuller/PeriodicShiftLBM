constexpr float resGlobal = 0.16f; 														// mm
constexpr int gridLevelCount = 3;
constexpr int wallRefinementSpan = 1;

constexpr int iterationCount = 150000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES
constexpr float SmagorinskyZoneLength = 10.f;

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number	
constexpr float rhoOutlet = 1.f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 16.94f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float RIn = 6.f;																// mm
constexpr float ROut = 12.f;															// mm
constexpr float C = 0.092784f;															// m2/s
constexpr int bladeCount = 5;

constexpr int iterationChunk = 2000;

#include "../include/types.h"

#include <sstream>

#include "../include/cellFunctions.h"
#include "../include/applyStreaming.h"
#include "../include/applyCollision.h"

#include "../include/boundaryConditions/applyBounceback.h"
#include "../include/boundaryConditions/applyMirror.h"
#include "../include/boundaryConditions/restoreRho.h"
#include "../include/boundaryConditions/restoreUxUyUz.h"
#include "../include/boundaryConditions/restoreRhoUxUyUz.h"
#include "../include/boundaryConditions/applyMBBC.h"

#include "../include/STLFunctions.h"
std::string STLPathBlade = "blade.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	float x, y, z;
	getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
	const float r2 = x * x + y * y;
	if ( r2 >= ROut * ROut )
	{
		Marker.bounceback = 1;
		return;
	}
	if ( r2 <= RIn * RIn )
	{
		Marker.bounceback = 1;
		return;
	}
	if ( kCell == 0 ) Marker.givenUxUyUz = 1;
	else if ( kCell == Info.cellCountZ-1 ) Marker.givenRho = 1;
	else Marker.fluid = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	float x, y, z;
	getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
	const float r = std::sqrt( x * x + y * y );
	const float vtPhys = (1.f / (r / 1000.f)) * C;
	const float vt = vtPhys * ( uzInlet / uzInletPhys );
	ux = - vt * (y / r);
	uy = vt * (x / r);
	uz = uzInlet;
	rho = rhoOutlet + Info.regulator;
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	const float zMax = Info.oz + Info.cellCountZ * Info.res;
	const float z = Info.oz + kCell * Info.res;
	if ( zMax - z < SmagorinskyZoneLength ) return 1.f; // 10mm from the end
	else return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const MarkerStruct &Marker, const InfoStruct &Info )
{
	rho = 1.f;
	ux = 0.f;
	uy = 0.f;
	uz = uzInlet;
	if (Marker.bounceback) uz = 0.f;
}

#include "../include/applyLocalCellUpdate.h"
#include "../include/plotter/exportSectionCutPlot.h"
#include "../include/fillEquilibrium.h"
#include "../include/gridRefinementFunctions.h"
#include "../include/DIADFunctions.h"
#include "../include/flowReportFunctions.h"

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

/*
float getVtDeviation( GridStruct &Grid, const float &r )
{
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	
	const int start = 0;
	const int end = (int)(3.14159f * 2.f * r / Info.res);
	const float evaluationZ = Info.oz + Info.cellCountZ * Info.res - SmagorinskyZoneLength;
		
	auto fetchVt = [ = ] __cuda_callable__( const int sinfgleIndex )
	{
		const float fi = ((float)sinfgleIndex / (float)end) * 2.f * 3.14159f;
		const float x = r * cosf(fi);
		const float y = r * sinf(fi);
		int iCell, jCell, kCell;
		getIJKCellIndexFromXYZ( iCell, jCell, kCell, x, y, evaluationZ, Info);
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		
		convertToPhysicalVelocity( ux, uy, uz, Info );
		
		const float vt = uy * cosf(fi) - ux * sinf(fi);
		
		return vt;
	};
	auto reductionVt = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	float vtSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchVt, reductionVt, 0.f );
	float vtAvg = vtSum / (float)end;
	const float vtTarget = 0.f; //(1.f / (r / 1000.f)) * C;
	const float vtDeviation = abs( vtTarget - vtAvg );
	return vtDeviation;
}
*/

/*
float getPLoss( GridStruct &Grid )
{
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	
	const int start = 0;
	const int end = Info.cellCountX * Info.cellCountY;
	
	auto fetchCount = [ = ] __cuda_callable__( const int sinfgleIndex )
	{
		const int iCell = sinfgleIndex % Info.cellCountX;
		const int jCell = sinfgleIndex / Info.cellCountX;
		const int kCell = 0;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		MarkerStruct Marker;
		Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0;
		else return 1;
	};
	auto reductionCount = [] __cuda_callable__( const int& a, const int& b )
	{
		return a + b;
	};
		
	auto fetchPLoss = [ = ] __cuda_callable__( const int sinfgleIndex )
	{
		const int iCell = sinfgleIndex % Info.cellCountX;
		const int jCell = sinfgleIndex / Info.cellCountX;
		const int kCell = 0;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		MarkerStruct Marker;
		Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		float pInlet = rho;
		convertToPhysicalPressure( pInlet, Info );
		
		float x, y, z;
		getXYZFromIJKCellIndex( iCell, jCell, kCell, x, y, z, Info );
		const float r = std::sqrt( x * x + y * y );
		const float vtPhys = (1.f / (r / 1000.f)) * C;
		
		const float pLoss = pInlet + 0.5f * rhoNominalPhys * vtPhys * vtPhys;
		
		return pLoss;
	};
	auto reductionPLoss = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	int count = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchCount, reductionCount, 0 );
	float pLossSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchPLoss, reductionPLoss, 0.f );
	float pLossAvg = pLossSum / (float)count;
	return pLossAvg;
}
*/

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

/*

void writeCaseResult( const std::vector<float>& historyVector, const int caseID ) {
	const int averagingIntervalStart = (iterationCount / 10 * 7);
	float sum = 0.f;
	int counter = 0;
	for ( int i = averagingIntervalStart; i < iterationCount; i++ )
	{
		sum += historyVector[i];
		counter++;
	}
	const float result = sum / (float)counter;
    std::ofstream outFile("optimizerResults.txt", std::ios::app);
    outFile << caseID << "; " << result << "\n";
    outFile.close();
}
*/

int main(int argc, char **argv)
{
	int caseID = 0;
    if (argc > 1) 
    {
		caseID = std::stoi(argv[1]);
	}

	STLStructCPU STLCPUBlade;
	readSTL( STLCPUBlade, STLPathBlade );
	STLStruct STLBlade( STLCPUBlade );
	checkSTLEdges( STLBlade );
	
	std::vector<STLStruct> STLs = { STLBlade };
	for ( int bladeIndex = 1; bladeIndex < bladeCount; bladeIndex++ )
	{	
		float radians = bladeIndex * 3.14159f * 2.f * (1.f / (float)bladeCount);
		STLStruct STLBladeRotated;
		STLBladeRotated = STLBlade;
		rotateSTLAlongZ( STLBladeRotated, radians );
		STLs.push_back( STLBladeRotated );
	}
	
	std::vector<DIADGridStruct> grids(gridLevelCount);
	// Coarse grid: Grid0
	grids[0].Info.res = resGlobal;
	grids[0].Info.dtPhys = dtPhysGlobal;
	grids[0].Info.nu = (grids[0].Info.dtPhys * nuPhys) / ((grids[0].Info.res/1000) * (grids[0].Info.res/1000));
	grids[0].Info.cellCountX = (int)( (2.f*ROut) / grids[0].Info.res ) + 2;
	grids[0].Info.cellCountY = grids[0].Info.cellCountX;
	grids[0].Info.cellCountZ = (STLBlade.zmax - STLBlade.zmin + SmagorinskyZoneLength + 20.f) / grids[0].Info.res;
	grids[0].Info.ox = - (grids[0].Info.cellCountX / 2) * grids[0].Info.res;
	grids[0].Info.oy = grids[0].Info.ox;
	grids[0].Info.oz = STLBlade.zmin - 10.f;
	grids[0].Info.cellCount = grids[0].Info.cellCountX * grids[0].Info.cellCountY * grids[0].Info.cellCountZ;
	buildIJKFromInfo( grids[0].IJK, grids[0].Info );
	
	buildDIADGrids( grids, STLs, 0 );
	
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
	
	std::vector<float> historyVector( iterationCount, 0.f );
	float uOutletMovingAvg = uzInlet;
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		updateGrid( grids, 0 );
		
		if (iteration%20 == 0 && iteration != 0)
		{
			FlowReportStruct FlowReport;
			XYZBoundsStruct Bounds;
			Bounds.xmin = -15.f;
			Bounds.xmax = 15.f;
			Bounds.ymin = -15.f;
			Bounds.ymax = 15.f;
			int kCut = 0;
			getFlowReportXY( grids, kCut, Bounds, FlowReport );
			float pAvg = FlowReport.rho;
			convertToPhysicalPressure( pAvg );	
			for ( int shifter = 0; shifter < 20; shifter++ )
			{
				historyVector[iteration-shifter] = pAvg;
			}
			
			// REGULATING OUTLET RHO TO STRONGLY DAMPEN ACUSTIC WAVES
			kCut = grids[gridLevelCount-1].Info.cellCountZ-1;
			getFlowReportXY( grids, kCut, Bounds, FlowReport );
			for ( int level = 0; level < gridLevelCount; level++ )
			{
				grids[level].Info.regulator = ( FlowReport.uz - uOutletMovingAvg ) / invSqrt3;
			}
			uOutletMovingAvg = 0.999 * uOutletMovingAvg + 0.001 * FlowReport.uz;
		}
		
		if ( iteration % iterationChunk == 0) 
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)cellUpdatesPerIteration * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			exportHistoryData( historyVector, iteration, 0 );
			
			int counter = 1;
			for (float r = RIn + 1.f; r < ROut; r = r + 1.f) 
			{
				exportSectionCutPlotToiletPaperZ( grids, r, iteration + 10*0 + counter );
				system("python3 ../include/plotter/plotterGridID.py");
				counter++;
			}	
			lapTimer.reset();
			lapTimer.start();	
		}
	}
	/*
	lapTimer.stop();
	auto lapTime = lapTimer.getRealTime();
	std::cout << "Finished iteration " << iterationCount << std::endl;
	const float updateCount = (float)Grid.Info.cellCount * (float)iterationCount;
	const float glups = updateCount / lapTime / 1000000000.f;
	std::cout << "GLUPS: " << glups << std::endl;
	
	int counter = 1;
	for (float r = RIn + 1.f; r < ROut; r = r + 1.f) 
	{
		exportSectionCutPlotToiletPaperZ( Grid, r, 10*0 + counter );
		system("python3 ../include/plotter/plotter.py");
		counter++;
	}
	
	exportHistoryData( historyVector, iterationCount, 0 );
	
	if (argc > 1) 
    {
		writeCaseResult( historyVector, caseID );
	}
	
	std::cout << "Finshed case " << caseID << std::endl;	
	*/
	return EXIT_SUCCESS;
}
