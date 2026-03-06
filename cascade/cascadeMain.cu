constexpr int caseID = 1;

constexpr float resGlobal = 0.1f; 														// mm
constexpr int iterationCount = 100000;
constexpr int iterationChunk = 10000;

constexpr float SmagorinskyConstantGlobal = 0.1f; 										// set to zero to turn off LES

constexpr float uzInlet = 0.05f; 														// also works as nominal LBM Mach number
constexpr float uInletAngle = 10.f;														// degrees
constexpr float uyInlet = (2.f * 3.14159f * uInletAngle / 360.f) * uzInlet;				// this gives hull angle		
constexpr float rhoOutlet = 1.f;
constexpr float nuPhys = 1e-6;															// m2/s water
constexpr float rhoNominalPhys = 1000.0f;												// kg/m3 water
constexpr float uzInletPhys = 20.f; 													// m/s
constexpr float dtPhysGlobal = (uzInlet / uzInletPhys) * (resGlobal/1000); 				// s

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

#include "../STLFunctions.h"
std::string STLPathVane = "vane.STL";

__cuda_callable__ void getMarkers( 	const int& iCell, const int& jCell, const int& kCell, 
									MarkerStruct &Marker, const InfoStruct& Info )
{
	if ( Marker.bounceback ) return;
	if ( kCell == 0 || kCell == Info.cellCountZ-1 ) Marker.periodic = 1;
	if ( jCell == 0 ) Marker.givenUxUyUz = 1;
	else if ( jCell == Info.cellCountY-1 ) Marker.givenRho = 1;
	else if ( iCell == 0 || iCell == Info.cellCountX-1 ) Marker.bounceback = 1;
	else Marker.fluid = 1;
}

__cuda_callable__ void getGivenRhoUxUyUz( 	const int& iCell, const int& jCell, const int& kCell, 
											float& rho, float& ux, float& uy, float& uz,
											InfoStruct& Info )
{
	ux = 0.f;
	uy = uyInlet;
	uz = uzInlet;
	rho = 1.f;
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
#include "../plotter/exportSectionCutPlot.h"
#include "../fillEquilibrium.h"

int main(int argc, char **argv)
{
	const float h = 10.f;
	const float a = 10.f;
	const float b = 3.f;
	const float l = 15.f;
	const float w = 1.f;
	const float t = 1.f;
	const float s = 15.f;
	
	constexpr float width = 34.f;
	
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
	Grid.Info.cellCountX = (int)( (width) / Grid.Info.res );
	Grid.Info.cellCountY = (int)( (5.f * h) / Grid.Info.res );
	Grid.Info.cellCountZ = (int)( (s) / Grid.Info.res );
	Grid.Info.ox = - 0.5f * Grid.Info.res * Grid.Info.cellCountX;
	Grid.Info.cellCount = Grid.Info.cellCountX * Grid.Info.cellCountY * Grid.Info.cellCountZ;
	Grid.fArray.setSizes( 27, Grid.Info.cellCount );
	Grid.shifter = IntArrayType( 27, 0 );
	fillEquilibriumFromFunction( Grid );
	
	Grid.bouncebackMarkerArray = BoolArrayType( Grid.Info.cellCount, 0 );
	const bool insideMarkerValue = 1;
	for ( int shift = -3; shift <= 3; shift++ )
	{
		STLVane.oz = 0.5f * s + shift * s;
		BoolArrayType bouncebackVane = BoolArrayType( Grid.Info.cellCount, 0 );
		applyMarkersInsideSTL( bouncebackVane, STLVane, insideMarkerValue, Grid.Info );
		sumBoolArrays( bouncebackVane, Grid.bouncebackMarkerArray, Grid.bouncebackMarkerArray );
	}
	std::cout << "Cell count: " << Grid.Info.cellCount << std::endl;
	
	//std::vector<float> historyMassFlow( iterationCount, 0.f );
	//std::vector<float> historyEta( iterationCount, 0.f );
	
	std::cout << "Starting simulation" << std::endl;
	
	TNL::Timer lapTimer;
	lapTimer.reset();
	lapTimer.start();
	for (int iteration=0; iteration<=iterationCount; iteration++)
	{
		applyStreaming( Grid );
		applyLocalCellUpdate( Grid );
		/*
		const int iStart = (-26.f - grids[gridLevelCount-1].Info.ox) / grids[gridLevelCount-1].Info.res;
		const int iEnd = (26.f - grids[gridLevelCount-1].Info.ox) / grids[gridLevelCount-1].Info.res;
		const int jStart = (-26.f - grids[gridLevelCount-1].Info.oy) / grids[gridLevelCount-1].Info.res;
		const int jEnd = grids[gridLevelCount-1].Info.cellCountY-1;
		
		float uzAvg, massFlow, pAvg;
		int iTemp, jTemp, flowReportK;
		const float xTemp = 0.f; 
		const float yTemp = 0.f;
		float flowReportZ = 30.f;
		getIJKCellIndexFromXYZ( iTemp, jTemp, flowReportK, xTemp, yTemp, flowReportZ, grids[gridLevelCount-1].Info);
		getFlowReport( flowReportK, grids[gridLevelCount-1], iStart, jStart, iEnd, jEnd, uzAvg, massFlow, pAvg );
		
		float uzAvgInlet, massFlowInlet, pAvgInlet;
		flowReportZ = -100.f;
		getIJKCellIndexFromXYZ( iTemp, jTemp, flowReportK, xTemp, yTemp, flowReportZ, grids[gridLevelCount-1].Info);
		getFlowReport( flowReportK, grids[gridLevelCount-1], 0, 0, grids[gridLevelCount-1].Info.cellCountX, grids[gridLevelCount-1].Info.cellCountY, uzAvgInlet, massFlowInlet, pAvgInlet );

		float lakePower = 0.5f * massFlow * uzAvgInlet * uzAvgInlet;
		float intakePower = 0.5f * massFlow * uzAvg * uzAvg + massFlow * (pAvg - pAvgInlet) / rhoNominalPhys;
		float eta = intakePower / lakePower;
		
		historyMassFlow[iteration] = massFlow;
		historyEta[iteration] = eta;
		*/
		
		if (iteration%iterationChunk == 0 && iteration!=0)
		{
			lapTimer.stop();
			auto lapTime = lapTimer.getRealTime();
			std::cout << "Finished iteration " << iteration << std::endl;
			const float updateCount = (float)Grid.Info.cellCount * (float)iterationChunk;
			const float glups = updateCount / lapTime / 1000000000.f;
			std::cout << "GLUPS: " << glups << std::endl;
			
			int iCut = Grid.Info.cellCountX / 2;
			exportSectionCutPlotZY( Grid, iCut, iteration );
			system("python3 ../plotter/plotter.py");
			
			//export3DPlot( grids[gridLevelCount-1], iteration + 30 );
			//system("python3 plotter3D.py");
			
			//exportHistoryData( historyMassFlow, historyEta, iteration, caseID );
			
			lapTimer.reset();
			lapTimer.start();
		}
	}
	std::cout << "Finshed successfuly" << std::endl;	
	return EXIT_SUCCESS;
}

/*
void getFlowReport( const int kCell, GridStruct &Grid, const int &iStart, const int &jStart, const int &iEnd, const int &jEnd,
						float &uzAvgPhys, float &massFlowPhys, float &pPhys )
{
	InfoStruct Info = Grid.Info;
	auto fArrayView  = Grid.fArray.getView();
	auto shifterView  = Grid.shifter.getConstView();
	bool useBouncebackArray = false;
	auto bouncebackMarkerArrayView = Grid.bouncebackMarkerArray.getConstView();
	if ( Grid.bouncebackMarkerArray.getSize() > 0 )
	{
		useBouncebackArray = true;
	}
	const int iSpan = iEnd - iStart;
	const int jSpan = jEnd - jStart;
	
	const int start = 0;
	const int end = iSpan * jSpan;
	
	auto fetchCellCount = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % iSpan + iStart;
		const int jCell = singleIndex / iSpan + jStart;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0;
		else return 1;
	};
	auto reductionCellCount = [] __cuda_callable__( const int& a, const int& b )
	{
		return a + b;
	};
	
	auto fetchUz = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % iSpan + iStart;
		const int jCell = singleIndex / iSpan + jStart;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0.f;
		
		int shiftedIndex[27];
		getShiftedIndex( cell, shiftedIndex, shifterView, Info );
		float f[27];
		for (int direction = 0; direction < 27; direction++) f[direction] = fArrayView( direction, shiftedIndex[direction] );	
		float rho, ux, uy, uz;
		getRhoUxUyUz(rho, ux, uy, uz, f);
		return uz;
	};
	auto reductionUz = [] __cuda_callable__( const float& a, const float& b )
	{
		return a + b;
	};
	
	auto fetchRho = [ = ] __cuda_callable__( const int singleIndex )
	{
		const int iCell = singleIndex % iSpan + iStart;
		const int jCell = singleIndex / iSpan + jStart;
		
		int cell;
		getCellIndex( cell, iCell, jCell, kCell, Info );
		MarkerStruct Marker;
		if ( useBouncebackArray ) Marker.bounceback = bouncebackMarkerArrayView( cell );
		getMarkers( iCell, jCell, kCell, Marker, Info );
		if ( Marker.bounceback ) return 0.f;
		
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
	
	const int cellCount = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchCellCount, reductionCellCount, 0 );
	float uzSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchUz, reductionUz, 0.f );
	float rhoSum = TNL::Algorithms::reduce<TNL::Devices::Cuda>( start, end, fetchRho, reductionRho, 0.f );
		
	float uzAvg = uzSum / (float)cellCount;
	float rhoAvg = rhoSum / (float)cellCount;
	float ux, uy = 0;
	convertToPhysicalVelocity( ux, uy, uzAvg, Info );
	convertToPhysicalVelocity( ux, uy, uzSum, Info );
	float pAvg = rhoAvg;
	convertToPhysicalPressure( pAvg );
	
	float volumetricFlow = uzSum * (Info.res / 1000.f) * (Info.res / 1000.f);
	float massFlow = volumetricFlow * rhoAvg * rhoNominalPhys;
	
	uzAvgPhys = uzAvg;
	massFlowPhys = massFlow;
	pPhys = pAvg;
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
*/
