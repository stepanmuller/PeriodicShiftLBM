constexpr float resGlobal = 0.3f; 														// mm
constexpr int iterationCount = 50000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES
constexpr float SmagorinskyZoneLength = 10.f;

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number	
constexpr float rhoOutlet = 1.f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 10.4f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

constexpr float RIn = 7.f;																// mm
constexpr float ROut = 15.f;															// mm
constexpr float C = 0.092784f;															// m2/s
constexpr int bladeCount = 5;

constexpr int iterationChunk = 25000;

#include "../types.h"

#include <sstream>

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
	rho = 1.f;
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
}

#include "../applyLocalCellUpdate.h"
#include "../plotter/exportSectionCutPlot.h"
#include "../fillEquilibrium.h"

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
	
	GridStruct Grid;
	Grid.Info.res = resGlobal;
	Grid.Info.dtPhys = dtPhysGlobal;
	Grid.Info.nu = (Grid.Info.dtPhys * nuPhys) / ((Grid.Info.res/1000) * (Grid.Info.res/1000));
	
	Grid.Info.cellCountX = (int)( (2.f*ROut) / Grid.Info.res ) + 2;
	Grid.Info.cellCountY = Grid.Info.cellCountX;
	Grid.Info.cellCountZ = (STLBlade.zmax - STLBlade.zmin + SmagorinskyZoneLength + 20.f) / Grid.Info.res;
	Grid.Info.cellCount = Grid.Info.cellCountX * Grid.Info.cellCountY * Grid.Info.cellCountZ;
	
	Grid.Info.ox = - (Grid.Info.cellCountX / 2) * Grid.Info.res;
	Grid.Info.oy = Grid.Info.ox;
	Grid.Info.oz = STLBlade.zmin - 10.f;
	
	Grid.fArray.setSizes( 27, Grid.Info.cellCount );
	Grid.shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( Grid );
	
	Grid.bouncebackMarkerArray = BoolArrayType( Grid.Info.cellCount, 0 );
	const bool insideMarkerValue = 1;
	for ( int bladeIndex = 0; bladeIndex < bladeCount; bladeIndex++ )
	{	
		float radians = 0.f;
		if ( bladeIndex > 0 ) radians = 3.14159f * 2.f * (1.f / (float)bladeCount);
		rotateSTLAlongZ( STLBlade, radians );
		BoolArrayType bouncebackBlade = BoolArrayType( Grid.Info.cellCount, 0 );
		applyMarkersInsideSTL( bouncebackBlade, STLBlade, insideMarkerValue, Grid.Info );
		sumBoolArrays( bouncebackBlade, Grid.bouncebackMarkerArray, Grid.bouncebackMarkerArray );
	}
	std::cout << "Cell count: " << Grid.Info.cellCount << std::endl;
	
	std::vector<float> historyVector( iterationCount+1, 0.f );
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		applyStreaming( Grid );
		applyLocalCellUpdate( Grid );
		
		int counter = 0;
		float vtDeviation = 0.f;
		for (float r = RIn + 0.5f; r < ROut; r = r + 0.5f) 
		{
			vtDeviation += getVtDeviation( Grid, r );
			counter++;
		}
		vtDeviation = vtDeviation / counter;
		float pLoss = getPLoss( Grid );
		
		//float F = log(1 + vtDeviation) + log(1 + std::max({0.f, pLoss}));
		float F = pLoss;
		historyVector[iteration] = F;
		
		if ( iteration > 0 && iteration % iterationChunk == 0) exportHistoryData( historyVector, iteration, caseID );
	}
	lapTimer.stop();
	auto lapTime = lapTimer.getRealTime();
	std::cout << "Finished iteration " << iterationCount << std::endl;
	const float updateCount = (float)Grid.Info.cellCount * (float)iterationCount;
	const float glups = updateCount / lapTime / 1000000000.f;
	std::cout << "GLUPS: " << glups << std::endl;
	
	int counter = 1;
	for (float r = RIn + 1.f; r < ROut; r = r + 1.f) 
	{
		exportSectionCutPlotToiletPaperZ( Grid, r, 10*caseID + counter );
		system("python3 ../plotter/plotter.py");
		counter++;
	}
	
	exportHistoryData( historyVector, iterationCount, caseID );
	
	if (argc > 1) 
    {
		writeCaseResult( historyVector, caseID );
	}
	
	std::cout << "Finshed case " << caseID << std::endl;	
	return EXIT_SUCCESS;
}
