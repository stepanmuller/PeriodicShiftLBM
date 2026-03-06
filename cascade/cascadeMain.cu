constexpr int caseID = 4;

constexpr float resGlobal = 0.05f; 														// mm
constexpr int iterationCount = 10000000;
constexpr int iterationChunk = 100000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float uxInlet = 0.05f; 														// also works as nominal LBM Mach number
constexpr float uInletAngle = 10.f;														// degrees
constexpr float uyInlet = (2.f * 3.14159f * uInletAngle / 360.f) * uxInlet;				// this gives incoming stream angle		
constexpr float rhoOutlet = 1.f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uxInletPhys = 20.f; 													// m/s
constexpr float dtPhysGlobal = (uxInlet / uxInletPhys) * (resGlobal/1000); 				// s

constexpr float invSqrt3 = 0.577350269f; 
constexpr float soundspeedPhys = invSqrt3 * (resGlobal/1000) / dtPhysGlobal; 			// m/s

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
#include "../boundaryConditions/periodicXY.h"

#include "../STLFunctions.h"
std::string STLPathVane = "vane.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	Marker.periodicX = 1;
	Marker.periodicZ = 1;
	if ( jCell == 0 ) Marker.givenUxUyUz = 1;
	else if ( jCell == Info.cellCountY-1 ) Marker.givenRho = 1;
	else Marker.fluid = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	ux = uxInlet;
	uy = uyInlet;
	uz = 0.f;
	rho = 1.f;
}

__cuda_callable__ float getSmagorinskyConstant( const int  &iCell, const int &jCell, const int &kCell, const InfoStruct &Info  )
{
	const float y = Info.oy + jCell * Info.res;
	if ( y > 30.f ) return 100.f;
	else return SmagorinskyConstantGlobal;
}

__cuda_callable__ void getInitialRhoUxUyUz( const int &iCell, const int &jCell, const int &kCell, float &rho, float &ux, float &uy, float &uz, const InfoStruct &Info )
{
	rho = rhoOutlet;
	ux = 0.f;
	uy = 0.f;
	uz = 0.f;
}

#include "../applyLocalCellUpdate.h"
#include "../plotter/exportSectionCutPlot.h"
#include "../fillEquilibrium.h"

float getAvgInletPressure( GridStruct &Grid )
{
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	
	const int start = 0;
	const int end = Info.cellCountX * Info.cellCountZ;
	
	auto fetchRho = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % Info.cellCountX;
		const int jCell = 0;
		const int kCell = singleIndex / Info.cellCountX;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		return rho;
	};
	auto reductionRho = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	float rhoSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchRho, reductionRho, 0.f );
	float rhoAvg = rhoSum / (float)end;
	float pAvg = rhoAvg;
	convertToPhysicalPressure( pAvg );
	return pAvg;
}

void exportHistoryData( const std::vector<float>& historyMassFlow, const std::vector<float>& historyEta, const int &currentIteration, int fileNumber ) {
    FILE* fp = fopen("/dev/shm/historyData.bin", "wb");
    if (!fp) return;
    
    int count = currentIteration + 1;
    fwrite(&count, sizeof(int), 1, fp);
    fwrite(historyMassFlow.data(), sizeof(float), count, fp);
    fwrite(historyEta.data(), sizeof(float), count, fp); // Append the ETA array
    fclose(fp);
    
    // Construct the command string to pass the number as an argument
    std::string cmd = "python3 historyPlotter.py " + std::to_string(fileNumber) + " &";
    system(cmd.c_str());
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

int main(int argc, char **argv)
{
	const float h = 10.f;
	const float a = 10.f;
	const float b = 2.f;
	const float l = 25.f;
	const float w = 1.f;
	const float t = 0.5f;
	const float s = 10.f;
	
	std::stringstream ss;
	ss << "python3 vaneGenerator.py " << h << " " << a << " " << b << " " << l << " " << w << " " << t << " " << s;
	std::system(ss.str().c_str());
	
	STLStructCPU STLCPUVane;
	readSTL( STLCPUVane, STLPathVane );
	STLStruct STLVane( STLCPUVane );
	checkSTLEdges( STLVane );
	STLVane.oy = h;
	
	GridStruct Grid;
	Grid.Info.res = resGlobal;
	Grid.Info.dtPhys = dtPhysGlobal;
	Grid.Info.nu = (Grid.Info.dtPhys * nuPhys) / ((Grid.Info.res/1000) * (Grid.Info.res/1000));
	
	std::cout << "Sizing domain" << std::endl;
	Grid.Info.cellCountX = (int)( (1*s) / Grid.Info.res );
	Grid.Info.cellCountY = (int)( (4.f * h) / Grid.Info.res );
	Grid.Info.cellCountZ = 1;
	Grid.Info.cellCount = Grid.Info.cellCountX * Grid.Info.cellCountY * Grid.Info.cellCountZ;
	Grid.fArray.setSizes( 27, Grid.Info.cellCount );
	Grid.shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( Grid );
	
	Grid.bouncebackMarkerArray = BoolArrayType( Grid.Info.cellCount, 0 );
	const bool insideMarkerValue = 1;
	for ( int shift = -5; shift <= 5; shift++ )
	{
		STLVane.ox = 0.5f * s + shift * s;
		BoolArrayType bouncebackVane = BoolArrayType( Grid.Info.cellCount, 0 );
		applyMarkersInsideSTL( bouncebackVane, STLVane, insideMarkerValue, Grid.Info );
		sumBoolArrays( bouncebackVane, Grid.bouncebackMarkerArray, Grid.bouncebackMarkerArray );
	}
	std::cout << "Cell count: " << Grid.Info.cellCount << std::endl;
	
	std::vector<float> historyVector( iterationCount, 0.f );
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		applyStreaming( Grid );
		applyBCPeriodicX( Grid );
		applyLocalCellUpdate( Grid );
		
		historyVector[iteration] = getAvgInletPressure( Grid );
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)Grid.Info.cellCount * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			int kCut = 0;
			const int plotNumber = 0;
			exportSectionCutPlotXY( Grid, kCut, plotNumber );
			system("python3 ../plotter/plotter.py");
			
			exportHistoryData( historyVector, iteration, caseID );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}
